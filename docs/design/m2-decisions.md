# M2 Design Decisions — Resolving the Open Questions

**Status:** design only. Nothing here is implemented.
**Companion to:** `docs/design/m2-sequence-ambiguity.md` (the main M2 design).
**Purpose:** record decisions for the five open questions in §11 of the main
design, plus two questions surfaced while writing the serology companion, so
that an M2 implementation can start from a settled contract. Each decision gives
the choice, the rationale, the rejected alternatives, and the consequences for
the engine and tests.

Decisions are labeled `D1…D8`. Where a decision changes or pins a field in the
main design, the affected section is cited.

---

## D1 — Alignment source of truth (main §11.1)

**Decision.** Build a **precomputed, versioned per-locus column model** as a
build-time artifact derived from the IPD-IMGT/HLA alignments
(`*_gen.msf` / `*_nuc.msf` and `hla.dat` feature coordinates). The engine reads
the precomputed artifact at runtime, never the raw MSF.

**Artifact contract.**

```
ColumnModel (one per locus per database_version):
    database_version: str                 # e.g. "IPD-IMGT/HLA 3.55.0"
    locus: str
    frame: "gen" | "nuc"                   # genomic vs cDNA alignment frame
    n_columns: int
    feature_of_column: tuple[FeatureKey, ...]   # length n_columns
    alleles: Mapping[Allele, PackedRow]    # per-allele base per column + status
    checksum: str                          # sha256 of the normalized artifact
```

`PackedRow` encodes, per column, one of `{A,C,G,T, gap, indeterminate}` plus the
feature status rollup (§4 of main design).

**Rationale.**

- **Reproducibility pinning.** The `checksum` + `database_version` go into
  `ObservationResult` metadata (QA §10.10). Two runs against the same pinned
  artifact are byte-identical regardless of IPD release churn or MSF parser
  changes. Parsing MSF live would couple results to the parser version.
- **Performance.** The column model is the hot input to signature extraction
  (main §8). Precomputing the gapped frame once avoids re-aligning per
  observation.
- **Auditability.** The artifact is a reviewable, diffable file; a scientist can
  inspect exactly which columns belong to which feature.

**Rejected.** Reading MSF directly at runtime — rejected for reproducibility and
startup-cost reasons. A live call to the IPD Feature Service — rejected: network
dependency, non-reproducible, unavailable offline.

**Consequences.**

- Add a `tools/build_column_model.py` (out of M2 engine scope; a data-prep
  utility) that emits one artifact per locus per pinned version.
- The engine depends on the artifact *interface*, not on IPD file formats. Tests
  use small hand-authored artifacts (a dozen synthetic alleles), so the test
  suite needs no large reference download.
- `reference_policy = strict` (main §4.0) checks are expressed against artifact status.

---

## D2 — Column model for indels / gaps (main §11.2)

**Decision.** Treat an alignment **gap as a first-class symbol distinct from a
base and distinct from indeterminate**, and make gap participate in signature
equality **only when the gap-bearing columns are inside the interrogated
region and have `PRESENT` status**.

Three-way meaning per column (extends the three-valued compare in main §4):

| Symbol | Meaning | In `compare_column()` (main §4.1) |
|--------|---------|----------------|
| `A/C/G/T` | a called base | normal COMPATIBLE / INCOMPATIBLE |
| `gap` (`-`) | an alignment gap that **is** sequenced as absent of base (a real deletion the assay observed) | `gap` vs `gap` = COMPATIBLE; `gap` vs base = INCOMPATIBLE |
| `indeterminate` (`?`) | position not observed / no reference | always INDETERMINATE (never confirms) |

The crucial distinction: a **sequenced deletion** (`gap`, PRESENT) is evidence
and can distinguish alleles; a **gap arising from missing data** must be encoded
as `indeterminate`, not `gap`, by the column-model builder (D1). This is a
builder invariant, asserted at artifact-build time.

**Rationale.** Indels are real, discriminating variation (e.g. a single-base
insertion creating a null). Collapsing gap into indeterminate would hide genuine
distinctions; collapsing indeterminate into gap would fake identity — the exact
error §4 forbids. Keeping them separate preserves both correctness directions.

**Rejected.** "Treat every gap as indeterminate" (loses real deletions);
"treat every gap as a fifth base unconditionally" (would turn missing-data gaps
into false evidence).

**Consequences.**

- `compare()` becomes a 3×-symbol table but stays O(1) per column.
- Builder (D1) must classify each gap as sequenced-deletion vs missing-data; the
  classification rule is itself pinned in the artifact provenance.
- New mandatory test: an indel-null allele (e.g. a known `*…:…N` caused by a
  frameshift insertion in exon 2) is **distinguishable** under an exon-2 assay
  and **indistinguishable** under an assay that does not read that column.

---

## D3 — Assay capability decomposition (supersedes `zygosity_aware`) (main §11.3, §2.3.1)

**Decision.** The single `zygosity_aware` set-vs-multiset flag is **removed** as
scientifically incorrect: for a clean diploid site the base *set* already
encodes zygosity (`|set| == 1` homozygous, `== 2` heterozygous), so a multiset
adds no information (main §2.3.1, §7.3). It is replaced by **four independent
capability properties**, each modeling a distinct physical phenomenon:

| Capability | Meaning | M2 status |
|-----------|---------|-----------|
| `detects_heterozygosity` | reports both bases at a het column | **true for M2** (observations complete) |
| `resolves_dosage` | counts copies (2:1 vs 1:1); meaningful only with >2 contributing chromosomes or CNV | inert in clean diploid M2; vocabulary reserved for the M3 mixture model |
| `allele_dropout` | stochastic loss of one allele's signal → false homozygosity | **not M2** → M3 |
| `resolves_phase` | cis/trans assignment across ≥2 het positions (replaces `phase_resolving`) | **modeled in M2**; default false |

**Rationale.** The old flag conflated (a) *information content of a complete
diploid observation* (deterministic, in M2) with (b) *stochastic dropout of the
second allele* (M3). Separating the four makes each testable and keeps M2
deterministic and complete. Only `resolves_phase` changes M2 output;
`resolves_dosage` is inert until mixtures/CNV (M3+); `allele_dropout` is M3.

**Rejected.** Keeping one boolean — rejected: it mislabels dropout as an
observation-model choice and wrongly claims a set cannot express zygosity.

**Consequences.**

- Sequence-assay config carries the four properties. In M2, validation requires
  `detects_heterozygosity = true` (an M2 assay observes completely) and
  `allele_dropout` absent/false (dropout is M3); `resolves_dosage` must be
  false/absent (no effect on clean diploid); `resolves_phase` is free.
- Named presets set only `resolves_phase` (SSO/short-read SBT false; long-read /
  cloned SBT true). No preset sets a dosage flag in M2.
- Test: the §7.3 example shows set == multiset for clean diploid; a separate M3
  test will exercise `resolves_dosage` under a 3-contributor mixture.

---

## D4 — Representation policy: exact by default (main §11.4, §3.3)

**Decision.** The default is **exact** representation; a superset is produced
**only on explicit configuration**. `representation_policy` has three levels:

| Policy | Behavior |
|--------|----------|
| `exact` (**default**) | Always emit the exact candidate set (`\|`-enumeration, or lossless `/` factorization when the set equals the Cartesian product). On a **grammatical limit** (main §3.3.1 — constraint the grammar cannot express) raise `RepresentationError`. On a **computational-expansion limit** (large but expressible), still emit the exact set — no silent truncation. |
| `superset` | On a grammatical limit or a configured expansion cap, emit the nearest safe **superset** and record a `RepresentationLimitation` whose `reason` distinguishes `grammar_gap` from `expansion_cap`. Never narrows. |
| `annotate_only` | Emit the superset in the GL String **and** attach the exact candidate set as structured metadata. |

**Rationale.** Review item 5: exact enumeration with `|` must be the default and
a conservative superset must require explicit opt-in. Grammatical limits (no
operator exists) are distinct from expansion limits (set is large); only the
former can justify a superset, and only under non-default policy. This replaces
the earlier "superset default."

In all cases: **never silently narrow**, and any approximation is recorded in
structured metadata (main §6.2). A superset is a conservative over-approximation
that contains the true candidate
set), so it can never exclude the truth — QA §10.5 is preserved under fallback.

**Rationale.** Different downstream uses have different tolerances. An imputation
benchmark may opt into a flagged superset; a method-validation study measuring
ambiguity precisely keeps the default `exact` so any grammatical limit surfaces
loudly as a `RepresentationError`; `annotate_only` serves tools that consume the
structured result rather than the GL text. One policy knob covers all three
without changing the engine core.

**Rejected.** "Superset by default" (the earlier draft) — rejected per review
item 5: it silently approximated when exact `|` enumeration was available.
"Always fail on any large set" — rejected: a large-but-expressible set is not a
correctness problem and must still be emitted exactly by default.

**Consequences.**

- `representation_policy` added to strategy schema, **default `exact`**.
- `RepresentationLimitation` records: the locus, the reason
  (`grammar_gap` | `expansion_cap`), the emitted superset, and
  (for `annotate_only`) the exact candidate pair set.
- Test: a grammatical-limit case raises under `exact`, emits a flagged superset
  under `superset`; a large-but-expressible case emits the exact `|` set under
  `exact` (no truncation). Assert superset ⊇ truth in the superset cases.

---

## D5 — Serology placement: M2 vs M2.5 (main §11.5)

**Decision.** Keep serology **in M2**, implemented as a **sibling observation
backend** that shares the pair-enumeration skeleton but swaps the signature
function and the reference data source. It is specified in the serology
companion (`docs/design/m2-serology.md`). Sequence and serology backends are
delivered together because they share the genotype-ambiguity algorithm (main
§3.2) and the output/metadata contract; splitting them would duplicate that
core or delay it.

**Rationale.** The thesis of M2 — ambiguity as a consequence of limited
observation, computed at the genotype level — is identical for both. Serology's
"observation" is a specificity multiset; sequence's is a base-column signature.
The *only* differences are the signature function, the reference artifact, and
the output vocabulary. Co-delivering keeps one algorithm, one test harness.

**Rejected.** Deferring serology to a later M2.5 — rejected: it would either
duplicate the pair-enumeration engine later or leave the engine's abstraction
unproven against a second backend (the best evidence the abstraction is right is
that a second, quite different backend plugs in cleanly).

**Consequences.**

- Introduce an `ObservationBackend` protocol (signature extraction + candidate
  reduction hooks) with two implementations: `SequenceBackend`,
  `SerologyBackend`. The genotype-ambiguity algorithm (main §3.2–§3.3) lives in
  the shared layer, not in either backend.
- Serology data (allele→specificity, broad/split lattice) is a separate pinned
  artifact, parallel to the column model (D1).

---

## D6 — Backend selection and the `ObservationBackend` boundary (new)

*Surfaced while writing D5; recorded for implementation clarity.*

**Decision.** `assay_class` selects the backend; the backend exposes exactly two
hooks to the shared genotype engine:

```
ObservationBackend (protocol):
    def signature(self, allele: Allele, config) -> Signature: ...
    def combined_observation(self, a1: Allele, a2: Allele, config) -> Observation: ...
    # equality of Observation is backend-defined (3-valued for sequence,
    # specificity-multiset for serology); the shared engine only calls
    # `combined_observation(x, y) == O`.
```

Routing:

| `assay_class` | backend |
|---------------|---------|
| `idealized` | M1 per-copy path (unchanged; not a backend) |
| `sso`, `sbt` | `SequenceBackend` |
| `serology` | `SerologyBackend` |

**Rationale.** The shared engine (pair enumeration, GL reduction, phase
regrouping, metadata) must not know whether it is comparing base columns or
serologic specificities. A two-method protocol is the minimal seam. It also
makes the §7.2 regression test backend-agnostic: the same assertion that `/`
cross-products are not emitted applies to serology.

**Consequences.**

- The shared engine is written once against `ObservationBackend`.
- Adding a future backend (e.g. an explicit-signature probe model, main §5
  `intermediate`) is a third implementation, no engine change.

---

## D7 — Null-allele semantics across backends (new)

*Surfaced reconciling main §10 (null-allele QA) with D2 and serology.*

**Decision.** A null allele's expression status is **not** a special case in the
engine; it is carried by whichever reference artifact the backend uses:

- **Sequence backend:** the null's causal variant is a column (a SNP, or an
  indel via D2). It is distinguishable iff that column is in the interrogated,
  PRESENT region — exactly the normal rule. A null caused by a variant **outside**
  the interrogated region is indistinguishable from its expressed sibling, which
  is the scientifically correct outcome.
- **Serology backend:** a null allele maps to a **blank/unresolved** specificity
  (it expresses no antigen). The combined observation of `expressed + null`
  therefore looks like a hemizygous/blank pattern, and the candidate set
  includes genotypes consistent with that blank — this is how serology's
  inability to see a null is modeled honestly.

**Rationale.** Treating nulls uniformly through the reference artifact (rather
than a hard-coded "nulls are special" branch) keeps the engine honest: whether a
null is observable is a property of the assay + reference, not an assertion the
engine makes because the truth name ends in `N`.

**Consequences.**

- Mandatory tests (both backends): a null distinguishable under one assay and
  indistinguishable under another; serology of `expressed + null` producing the
  blank-consistent candidate set.
- Reconciles with M1's rule that an expression suffix is never migrated across
  truncation (`Allele.truncate_to_fields`): M2 never *infers* a null from the
  name; it reads expression from reference data.

---

## D8 — Phase-group enumeration is a product of molecules, not of loci (main §6.4, §6.4.1)

**Decision.** For loci in a declared phase group backed by a **phase-resolving**
assay, candidate genotypes are enumerated as **haplotype pairings over the
joint interrogated signature** of the grouped loci — i.e. the product of the
observed *molecules* — and never as the Cartesian product of the independent
per-locus genotype ambiguities with phase re-attached. Per-locus enumeration
followed by unrestricted regrouping is **forbidden** for phase-resolved groups.

**Rationale.** Review item 6. Even when each per-locus ambiguity is individually
correct, their unrestricted combination invents cross-locus haplotypes the assay
never observed (worked, verified, in main §6.4.1: 8 naive vs 4 correct
diplotypes; the fabricated `B*A~C*B + B*2~C*1` is excluded by the joint
observation). The molecule is the unit the phase-resolving assay actually reads;
the candidate set must be built from molecules.

**Two regimes, made explicit:**

| Assay phase capability over the group | Enumeration | GL output |
|---------------------------------------|-------------|-----------|
| `resolves_phase = false` (no joint observation) | independent per-locus; **no** regroup into one block | separate `^` groups; asserts no interlocus pairing (main §6.4 Rule A) |
| `resolves_phase = true` (joint molecule observation) | product of observed molecules over the joint signature | one `LocusGroup` spanning the grouped loci; `~` haplotypes, `\|` only between molecule-consistent pairings (main §6.4.1) |

**Rejected.** "Enumerate per locus, then take the Cartesian product and attach
phase from a declared phase group" — rejected: it fabricates haplotypes (main
§6.4.1) and would also silently reintroduce truth phase that was not observed
(violates spec §6 / QA §10.2).

**Consequences.**

- The engine has two code paths: independent per-locus (default), and
  joint-signature enumeration for phase-resolved groups. The joint path applies
  the §3.2 pair-enumeration algorithm to the *concatenated* multi-locus
  signature of the grouped loci.
- M1's `_assemble_groups` regrouping is retained **only** for the independent
  path; it must never be used to manufacture a phased block from independently
  enumerated per-locus candidates.
- Tests: `test_joint_phase_product_of_molecules_not_loci`,
  `test_unphased_regroup_asserts_no_interlocus_pairing`,
  `test_phase_regroup_never_invents_haplotype` (main §10).

**CONFIRMED (approved final decision).** The whole-molecule phase model is
approved **as an idealized assay model**. Two binding constraints follow:

1. **Joint observation must be declared explicitly.** The strategy specifies
   *which loci and which interrogated regions are jointly observed on the same
   molecule* — a `joint_observation_group` naming the member loci and their
   shared interrogated region. The engine enumerates molecule pairings only over
   a declared group.
2. **`resolves_phase = true` does not imply global phase.** It is a per-assay
   capability meaning "this assay *can* resolve phase within a declared joint
   group," **not** "phase is established across every typed locus." A locus with
   `resolves_phase = true` but no `joint_observation_group` membership is still
   reported with no interlocus phase. Phase is established only for loci that are
   (a) in a declared joint-observation group and (b) truly phased in truth.

This prevents the engine from silently asserting genome-wide phase just because
an assay is phase-capable. It also replaces any notion that a declared M1
`phase_group` alone (without a joint-observation region) triggers molecule
enumeration — a `phase_group` expresses *intent to retain* phase; a
`joint_observation_group` expresses *what the molecule physically spans*. Both
must be present to enumerate molecule pairings.

---

## Settled contract summary (for the implementer)

- Reference input = pinned, checksummed **column-model artifact** (D1); serology
  = parallel pinned specificity artifact (D5/D7).
- Gaps vs indeterminate are **distinct symbols** (D2).
- `zygosity_aware` is **removed**; replaced by four capability properties
  (`detects_heterozygosity`, `resolves_dosage`, `allele_dropout`,
  `resolves_phase`); only `resolves_phase` affects M2 (D3).
- `representation_policy ∈ {exact(default), superset, annotate_only}` (D4);
  exact `|` enumeration by default, superset only on opt-in, never narrowing.
- `reference_policy ∈ {strict, permissive(default), sequence_supplied}` (main
  §4.0); permissive results labeled `APPROXIMATE_INDETERMINATE`; incomplete-
  reference alleles stay in the candidate universe (main §4.2).
- Candidate comparison is **three-valued** (CONFIRMED_COMPATIBLE / INCOMPATIBLE
  / INDETERMINATE), not Boolean (main §4.1).
- Serology ships **in M2** as a second `ObservationBackend` (D5/D6).
- Backend seam = two methods: `signature`, `combined_observation` (D6).
- Nulls flow through reference data, no special engine branch (D7).
- Phase-resolved groups enumerate the **product of observed molecules**, not of
  loci; independent per-locus regrouping is forbidden for them (D8).

None of this is implemented; it is the agreed basis for the M2 implementation
once review signs off.
