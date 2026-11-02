# Milestone 2 Design — Sequence-Derived Ambiguity Engine

**Status:** design only. Nothing in this document is implemented yet.
**Audience:** TissueTypeForge maintainers (SIP / HLA Forge).
**Scope:** spec §4 (scientifically realistic ambiguity), §5 (sequence-based
second-field, intermediate, full-gene), §8 (GL String fidelity), the relevant
QA items in §10, and the M2 milestone in §11.

This document specifies *how* a configured sequence-based assay turns known
truth into the genotype it would actually report, where the report is limited
by **which sequence the assay interrogates** and by the fact that a **diploid
assay observes two chromosomes at once**. It deliberately does **not** reuse the
legacy "look the true allele up in a MAC table" shortcut (spec §4).

It builds directly on the Milestone 1 AST and contracts (see
`src/tissuetypeforge/glstring/model.py`, `observation/framework.py`,
`result.py`). The reader is assumed familiar with those.

---

## 0. Design goals and non-goals

**Goals**

1. Model ambiguity as a *consequence of limited observation*, derived from
   reference sequence, not as a lookup.
2. Treat the diploid observation as a single combined signature over two
   chromosomes, and compute **genotype-level** ambiguity from it (spec §4).
3. Keep allele ambiguity (`/`) and genotype ambiguity (`|`) strictly separate
   and only collapse one into the other when it is provably valid.
4. Handle missing / incomplete reference sequence explicitly — never treat
   "no data" as "sequence matches" (spec §4).
5. Emit a GL String that is exactly as precise as the observation justifies,
   and report explicitly when the representation cannot express the result
   faithfully (spec §8).

**Non-goals (deferred to M3+)**

- Stochastic error (dropout, false homozygosity, substitution): M3.
- Serology mappings and intermediate *probe* assays are specified structurally
  here but their reference data loading is a sibling concern (§5 of spec); this
  document focuses on the sequence engine. Serology is sketched in §9.
- Registry-scale batching: M4.

---

## 1. Where the engine sits in the pipeline

M1 pipeline (unchanged):

```
truth GL String -> parse -> Diplotype(truth)
  -> per locus: coverage + resolution -> observed per-copy AlleleLists
  -> phase regrouping -> observed Diplotype -> serialize
```

M1 resolves each **copy independently** (`observation/resolution.py`
`resolve_allele_list` is applied per copy). That is correct *only* for the
idealized modes, where the method is asserted to resolve each copy exactly.

**M2 replaces the per-locus step for sequence-based assays with a genotype-level
observation.** The unit of observation becomes the *unordered pair of true
alleles at a locus* (plus their phase identity, carried through unchanged), not
each allele separately. The output of the engine for a locus is a
`LocusGroup` whose `alternatives` are the candidate `GenotypeBlock`s — i.e.
genotype ambiguity (`|`) is produced here. Phase regrouping and serialization
downstream are reused as-is.

```
truth Diplotype(locus) = unordered {A_true1, A_true2}
  -> signature extraction (per assay)       §2
  -> observed combined signature            §2.3
  -> candidate allele set per chromosome    §3.1
  -> genotype reconstruction (pairing)      §3.2   <-- the core of M2
  -> LocusGroup with | alternatives         §3.3 / §6
```

---

## 2. Observable sequence model

### 2.1 Reference substrate

An allele's observable content is derived from a **versioned reference
sequence**, not from its name. The engine needs, per `(database_version, allele)`:

- the genomic/cDNA sequence partitioned into **named features** (exons,
  introns, UTRs) with coordinates, following the IPD-IMGT/HLA feature model
  (see `hla.dat` / the Feature Service / `hla_nom` + `_gen.msf`/`_nuc.msf`
  alignments);
- a flag per feature indicating whether sequence for that feature is **present,
  partial, or absent** for that allele (critical for §4).

This is captured by a `ReferenceSequence` abstraction (data loading is out of
scope for this design but the interface is fixed here):

```
AlleleSequence:
    allele: Allele
    database_version: str
    features: Mapping[FeatureKey, FeatureSequence]

FeatureKey:            e.g. ("exon", 2), ("intron", 1), ("utr", 3)
FeatureSequence:
    status: PRESENT | PARTIAL | ABSENT
    sequence: str | None         # None when ABSENT
    span: (start, end) | None    # coordinates in the alignment frame
```

All sequences for one locus must be expressed in a **common alignment frame**
(gapped alignment columns) so that two alleles can be compared position-by-
position. The engine operates on alignment columns, never on raw unaligned
indices.

### 2.2 Interrogated region = what the assay can see

A sequence assay is configured (extending the M1 `LocusConfig`) with:

- `interrogated_features`: the set of `FeatureKey`s the assay reads
  (e.g. exon-2+exon-3 class I; exon-2 class II; all features for full-gene).
- `interrogation_mask` (optional): a finer per-column mask within a feature,
  for assays that do not read a whole feature cleanly (e.g. primer-limited
  amplicons). Defaults to the whole feature.
- `resolves_phase`: whether the assay resolves cis/trans *within the
  interrogated region* (almost always false for short-read / SSO; may be true
  for long-read SBT). One of the four capability properties of §2.3.1; feeds
  the genotype-pairing step, §3.2, and the joint-phase evaluation, §6.4.

The **observable signature** of an allele under an assay is the projection of
that allele's reference sequence onto the interrogated columns:

```
signature(allele, assay) = tuple of (column -> base) over interrogated columns
                           for features whose status is PRESENT (see §4 for
                           PARTIAL / ABSENT handling)
```

Two alleles are **indistinguishable under the assay** iff their signatures are
equal *and* neither has missing data in the interrogated region (precise
definition in §4).

### 2.3 Diploid combined observation — the central modeling choice

A diploid assay does not read two separate haplotype sequences. It reads a
**mixture**: at every interrogated column it observes the bases present across
both chromosomes.

```
observed_column(A1, A2, c) = { base(A1, c), base(A2, c) }    # a SET of bases
combined_observation(A1, A2, assay) = tuple of observed_column over interrogated c
```

Two important properties follow, and they are the reason M2 cannot match
alleles independently:

1. **Heterozygous columns lose chromosome assignment.** Observing `{G, T}` at a
   column tells you one chromosome has G and the other T, but not *which*. When
   several columns are heterozygous, the assignment of variants across columns
   to chromosomes is ambiguous (classic *phase ambiguity* of unphased SBT).
2. **The observation is a property of the pair, not of each allele.** The same
   column-wise base sets can be produced by more than one *pair* of alleles.
   That is genotype ambiguity and must be enumerated over pairs.

#### 2.3.1 Zygosity is already in the set — do not conflate four distinct capabilities

An earlier draft carried a single `zygosity_aware` set-vs-multiset flag. That
was **scientifically wrong** and is removed. For a **clean diploid** site with
exactly two contributing chromosomes, the *set* of observed bases already
encodes zygosity:

```
|observed_column| == 1   <=>   homozygous at that column   (e.g. {G})
|observed_column| == 2   <=>   heterozygous at that column  (e.g. {G,T})
```

A diploid site can show at most two distinct bases, so a **multiset** (dosage)
adds **no information** over a set at a single clean-diploid column: there is no
third copy to count. The set already distinguishes `{G}` (homozygous G) from
`{G,T}` (heterozygous). Asserting that "set semantics cannot tell homozygous
from heterozygous" confused *a qualitative assay's technical dropout of the
second allele* (an M3 stochastic phenomenon) with *the information content of a
complete diploid observation* (deterministic, M2).

We therefore decompose what was one flag into **four independent capability
properties** on the assay config, each with a precise, separate meaning:

| Capability | Meaning | M2 (deterministic/complete) | Deferred |
|-----------|---------|------------------------------|----------|
| `detects_heterozygosity` | the assay reports both bases at a het column (i.e. observes the full set) | **Assumed true in M2** for interrogated, PRESENT columns — M2 observations are complete | — |
| `resolves_dosage` | the assay counts copies (2:1 vs 1:1), relevant only when >2 chromosomes contribute (mixture) or under copy-number variation | **Out of scope in clean diploid M2**; a multiset/count model is only meaningful for mixtures | mixture/contamination → M3 |
| `allele_dropout` | a true allele's signal is stochastically lost, collapsing an apparent het to a false hom | **Not in M2** (M2 is complete and deterministic) | → M3 |
| `resolves_phase` | cis/trans assignment of variants across ≥2 heterozygous positions is determined (long-read / cloned SBT) | **Modeled in M2** via `resolves_phase` (replaces the old `phase_resolving`); default false | — |

Consequences for M2:

- **M2 observations are complete and deterministic.** Every interrogated,
  PRESENT column contributes its full base set. Heterozygosity is therefore
  always detected at observed columns; there is no "qualitative assay hides the
  second allele" case in M2 — that is precisely allele dropout and lives in M3.
- **Dosage/`resolves_dosage` is not exercised by clean diploid truth.** It is
  retained in the config vocabulary only so the M3 mixture model and any future
  CNV handling (e.g. DRB copy number) can switch to a multiset observation
  without a schema change. In M2 it must be left at its default and does not
  alter the candidate set.
- **`resolves_phase`** is the only one of the four that changes M2 output; it
  feeds the intralocus cis/trans decision in §3.2 and §7.2 and the joint
  evaluation in §6.4.

See §7.3 for an analytically valid example that states its assumptions
explicitly and shows the *only* regime in which a multiset would differ from a
set (a >2-chromosome mixture), confirming it is out of M2 scope.

---

## 3. Diploid genotype ambiguity algorithm

Input for one locus: the true unordered pair `{A1, A2}` (from truth), the assay
config, and the reference database. Output: a `LocusGroup` (one or more
candidate `GenotypeBlock`s joined by `|`).

### 3.1 Candidate universe

Let `U` be the set of alleles of this locus in the configured database that have
**usable** sequence over the interrogated region (status rules in §4). For each
`u ∈ U` compute `signature(u, assay)`.

### 3.2 Genotype reconstruction (pair enumeration) — do NOT match independently

Compute the *true* combined observation once:

```
O = combined_observation(A1, A2, assay)
```

Then enumerate the **unordered pairs** `{x, y}` with `x, y ∈ U` such that the
pair *reproduces the same observation* under the three-valued comparison of §4:

```
candidates = { {x, y} : compare( combined_observation(x, y, assay), O ) ≠ INCOMPATIBLE }
```

(The comparison is three-valued — CONFIRMED_COMPATIBLE / INCOMPATIBLE /
INDETERMINATE — not Boolean equality; see §4. A pair is a candidate unless it is
positively INCOMPATIBLE; INDETERMINATE pairs are kept and flagged.)

This is the definition of genotype ambiguity: every pair not excluded by what
the assay observed. Properties:

- `{A1, A2}` is always in `candidates` (truth reproduces its own observation) —
  this guarantees QA §10.5 (ambiguity set contains the truth when no error).
- A pair `{x, y}` can be present even if **neither `x` nor `y` individually
  equals `A1` or `A2`**. This is exactly the case independent matching gets
  wrong (worked in §7.1–§7.2).
- Homozygous truth (`A1 == A2`) is handled naturally: every observed column is
  a singleton set `{b}`, which positively **excludes** any heterozygous pair
  (whose column would be a 2-element set) — because M2 observations are complete
  (§2.3.1), a singleton column is a confirmed homozygous observation, not an
  "allele hidden" one. Candidate homozygotes are other alleles whose
  interrogated sequence is identical to truth's; candidate heterozygotes arise
  only where the reference is INDETERMINATE (§4, §7.4), never from set-vs-
  multiset ambiguity.

Complexity: naive enumeration is `O(|U|^2)` per locus. This is acceptable for a
reference design; §8 lists the optimization (signature bucketing) that reduces
it to grouping by signature and pairing buckets, i.e. roughly
`O(|U| + k^2)` where `k` is the number of *distinct* signatures involved.

### 3.3 Reduction to a minimal GL String form

`candidates` is a set of unordered allele pairs. **The default is exact
representation.** The engine serializes the exact candidate set and only ever
produces a superset when explicitly configured to (review item 5). Reduction
lattice, in order:

1. **Single unambiguous genotype** — if `|candidates| == 1`. Serialized as one
   `GenotypeBlock` (`x+y`).
2. **Allele-ambiguity factorization (`/`)** — collapse to
   `(/-list of L) + (/-list of R)` **only if** `candidates` equals the full
   unordered Cartesian product `L × R` (every combination present). This is a
   purely cosmetic, lossless rewrite of an exact set. If the set does not
   factor exactly, **do not** force `/`.
3. **Exact genotype ambiguity (`|`)** — otherwise, emit the exact irreducible
   set of `GenotypeBlock` alternatives joined by `|`, with canonical ordering
   (sort alleles within a pair, sort pairs) for deterministic, round-trippable
   output. **This is the default for any non-factorable set — prefer exact `|`
   enumeration over any approximation.**

#### 3.3.1 Grammatical limits vs computational expansion limits

These are different and must not be conflated:

- **Grammatical limit:** the GL String grammar *cannot express* the constraint
  at all (e.g. a cross-locus joint constraint, or a partial intralocus phase
  that is neither full cis/trans nor full ambiguity — §6.3). Exact `|`
  enumeration does not help because the needed relation has no operator. Only
  here is a flagged superset ever appropriate.
- **Computational expansion limit:** the exact `|` set is expressible but
  *large* (many alternatives). This is a performance/size concern, **not** a
  correctness one. The default is still to emit it exactly. Truncating it to a
  superset for size reasons requires explicit opt-in.

#### 3.3.2 Superset only on explicit configuration

A conservative superset (a `|` set, or `/`-cross-product, that strictly contains
the true candidate set) is emitted **only** when the strategy sets
`representation_policy` to allow it (decision D4):

- `exact` (**default**): never emit a superset. On a grammatical limit, raise
  `RepresentationError`. On a computational-expansion limit, still emit the
  exact set (no silent truncation).
- `superset`: on a grammatical limit *or* a configured expansion cap, emit the
  nearest safe superset and record a `RepresentationLimitation` with the
  distinguishing `reason` (`grammar_gap` vs `expansion_cap`).
- `annotate_only`: emit the superset in the GL String but attach the exact
  candidate set as structured metadata.

In every case: **never silently narrow**, and record any approximation in
structured metadata (§6.2). The default path produces an exact, non-approximate
GL String.

#### 3.3.3 Configurable resource limits produce explicit errors

Approved final decision on D4: exact enumeration is default, but it must be
**bounded**. The engine exposes configurable resource limits; exceeding one
raises an **explicit error**, never unbounded enumeration and never a silent
approximation.

| Limit | Meaning | On breach |
|-------|---------|-----------|
| `max_candidate_pairs` | cap on enumerated candidate genotype pairs per locus/group | raise `ResourceLimitError(kind="candidate_pairs")` |
| `max_alternatives` | cap on `\|` alternatives in the emitted GL String | raise `ResourceLimitError(kind="alternatives")` |
| `max_universe` | cap on candidate-universe size `|U|` loaded for a locus | raise `ResourceLimitError(kind="universe")` |
| `max_joint_molecules` | cap on molecule pairings enumerated for a joint-phase group (§6.4.1) | raise `ResourceLimitError(kind="joint_molecules")` |

Interaction with `representation_policy`:

- Under `exact` (default), a resource breach is a **hard error** — the engine
  will not substitute a superset to stay within budget (that would be a silent
  approximation).
- Under `superset`/`annotate_only`, the strategy *may* additionally set
  `on_resource_limit = "superset"` to convert an `expansion_cap`/`alternatives`
  breach into a flagged superset (reason `expansion_cap`) instead of an error.
  This conversion is explicit and recorded; it is never the default.

The SIP `glstring` library already raises `ExpansionLimitError` for its own
canonicalization/expansion operations; the engine maps that onto
`ResourceLimitError` so there is one limit-breach contract surfaced to callers.

- The engine produces, **per locus**, a `LocusGroup` of candidate pairs over
  that single locus (§3.2).
- **Interlocus phase regrouping (`_assemble_groups`, M1) runs afterwards and
  must not invent cross-locus allele associations.** When two loci are NOT
  jointly phase-observed, their per-locus candidate sets are reported as
  independent groups; the engine never forms the cross-product of one locus's
  candidates with another's as if they were linked. See §6.4 for the formal
  rule and a synthetic example.
- Phase is still only *retained* when truth phased the loci **and** a phase
  group is declared. M2 adds no new interlocus phase.
- **Intralocus** phase (cis/trans of variants within one locus) uses the
  `resolves_phase` capability (§2.3.1, renamed from `phase_resolving`): if
  `resolves_phase` is false and the interrogated region has ≥2 heterozygous
  columns whose cross-assignment yields distinct allele pairs, those pairs
  appear as separate `|` alternatives (§7.2). If `resolves_phase` is true, the
  cross-assignment is fixed and those alternatives are excluded.

---

## 4. Three information layers; missing / incomplete reference

Correction (review item 2): the engine must keep three distinct layers separate
and never collapse them.

1. **True biological sequence** of the individual's two alleles. This is what
   the laboratory assay physically interrogates. It exists whether or not IPD
   has it on record.
2. **Available IPD-IMGT/HLA reference sequence** for each named allele. This is
   our *knowledge*, and it may be incomplete (PARTIAL / ABSENT features).
3. **The assay's observable information** = the true biological sequence
   projected onto the interrogated region. This is what an assay *reports*.

The critical rule: **incompleteness of the reference database (layer 2) is not
absence of signal in the laboratory (layer 3).** When IPD lacks a base for a
*candidate* allele, that is a gap in our ability to compare the candidate to the
observation — it must become INDETERMINATE, never a silent MATCH and never a
silent exclusion.

The two layers where "unknown" can arise are therefore different in kind:

- **Observed side (truth projected):** a column the *assay* did not read is
  simply not in the interrogated region — it contributes nothing. A column the
  assay read is always known for the true alleles (M2 is complete, §2.3.1).
- **Candidate side (reference knowledge):** a candidate allele may have no IPD
  base at an interrogated column. That is the INDETERMINATE case below.

For a given interrogated column `c`, comparing the observed base(s) against a
candidate allele `u`:

- `u` has `PRESENT` reference at `c`: candidate base is known.
- `u` has `PARTIAL`/`ABSENT` reference covering `c`: candidate base is **unknown
  from reference** (`?`) — a layer-2 gap, not a layer-3 absence.

Comparison is three-valued (see §4.1 for the formal result type):

```
compare_column(observed_set, candidate_base):
    if candidate_base == '?':   -> INDETERMINATE   # reference gap (layer 2)
    elif candidate_base ∈ observed_set:  -> COMPATIBLE
    else:                       -> INCOMPATIBLE     # positive contradiction
```

Rules:

1. **An INCOMPATIBLE column excludes the pair.** A candidate base positively
   absent from the observed set at a read column is hard evidence.
2. **An INDETERMINATE column cannot confirm compatibility.** A pair compatible
   on all PRESENT columns but INDETERMINATE elsewhere is kept as an
   *indeterminate candidate*, flagged `uses_incomplete_reference = True`; the
   affected alleles are recorded. The engine never promotes INDETERMINATE to
   COMPATIBLE.
3. **If the true allele's reference is itself incomplete over the interrogated
   region**, fewer columns discriminate, generally producing a *larger*
   candidate set — the honest outcome.

### 4.0 strict / permissive / sequence-supplied behavior

The assay config selects how layer-2 gaps are handled via a required
`reference_policy`. **Default is `strict`** (approved final decision).

| Policy | Behavior when an INDETERMINATE (layer-2 reference gap) column affects whether a candidate belongs in the exact set | Result labeling |
|--------|-------------------------------------------------------------------------------------------------------------------|-----------------|
| `strict` (**default**) | The engine returns the exact candidate set **only when it can prove exactness**. If any candidate's membership is INDETERMINATE because of incomplete reference (its inclusion/exclusion cannot be decided), `strict` **raises `IndeterminateReferenceError`** naming the affected alleles. It must **never** silently drop indeterminate candidates and present the remainder as exact. Also raises if the **true** allele lacks reference over the interrogated region. | exact — or an explicit error; never a silently narrowed "exact" set |
| `permissive` | Keep INDETERMINATE candidates in the set (§4.1 rule 2). The result is explicitly an **approximate / indeterminate** set, never an exact observation. | **labeled `APPROXIMATE_INDETERMINATE`** in metadata; indeterminate candidates tagged individually; GL String is a superset of the complete-reference result |
| `sequence_supplied` | The caller supplies observed sequence for the interrogated region (truth side), so truth-side gaps vanish. **Candidate-side uncertainty from incomplete reference does not vanish** — those candidates are still INDETERMINATE and are governed by the companion policy (`strict` → error, `permissive` → flagged). | exact for the supplied truth region; candidate-side gaps as per companion policy |

Key requirements (approved final decision on reference policy):

1. **`strict` is the default** and is *fail-closed*: when incomplete reference
   knowledge prevents determining the exact candidate genotype set, it produces
   an **explicit error** (`IndeterminateReferenceError`), not a quietly reduced
   set.
2. **No silent exclusion-and-call-exact.** The forbidden behavior of the earlier
   draft (drop indeterminate candidates, label the rest exact) is removed.
3. **`permissive` is an explicit opt-in** producing flagged approximate results
   (`APPROXIMATE_INDETERMINATE`).
4. **`sequence_supplied` does not eliminate candidate-side uncertainty.**
   Supplying truth sequence removes truth-side gaps only; incomplete reference
   for *candidate* alleles still yields INDETERMINATE handled by the companion
   policy.

### 4.1 Three-valued candidate comparison (replaces Boolean equality)

Correction (review item 4): the engine must **not** test `combined_observation(x,y)
== O` with Boolean equality. It computes a three-valued result per candidate
pair, aggregated over interrogated columns:

```
CandidateComparison = CONFIRMED_COMPATIBLE | INCOMPATIBLE | INDETERMINATE

compare_pair(O, candidate_pair):
    results = [ compare_column(O[c], candidate base set at c) for c in columns ]
    if any result is INCOMPATIBLE:            -> INCOMPATIBLE
    elif any result is INDETERMINATE:         -> INDETERMINATE
    else:                                     -> CONFIRMED_COMPATIBLE
```

(Where a candidate pair's per-column contribution is itself a set
`{base(x,c), base(y,c)}` compared to the observed set `O[c]`: the pair is
column-INCOMPATIBLE if the two sets differ with no reference gap, column-
INDETERMINATE if either side has a `?`, else column-COMPATIBLE.)

**Effect on the returned ambiguity set:**

- `CONFIRMED_COMPATIBLE` pairs form the **core** candidate set. These are
  reported as genotype ambiguity regardless of policy.
- `INCOMPATIBLE` pairs are excluded (hard evidence).
- `INDETERMINATE` pairs are handled by `reference_policy` (§4.0): included and
  flagged under `permissive`/`sequence_supplied`, excluded under `strict`.

The returned `LocusGroup` therefore has a precise provenance: each `|`
alternative is tagged `confirmed` or `indeterminate`. Metadata (§6.2) reports
`confirmed_candidate_count` and `indeterminate_candidate_count` separately so a
downstream consumer can distinguish "the assay genuinely cannot resolve these"
(confirmed ambiguity) from "our reference knowledge is too incomplete to decide"
(indeterminate ambiguity). The truth pair is always at least INDETERMINATE and,
when truth's own reference is PRESENT over the interrogated region,
CONFIRMED_COMPATIBLE — preserving QA §10.5.

### 4.2 Incomplete-reference alleles stay in the candidate universe

The candidate universe `U` (§3.1) is **all alleles of the locus in the pinned
database**, including those with PARTIAL/ABSENT reference over the interrogated
region. An allele is **never silently dropped from `U` merely because its
reference is incomplete.** Its incompleteness surfaces as INDETERMINATE
comparisons (§4.1), governed by policy (§4.0), and is recorded in metadata.
Silently excluding reference-incomplete alleles would fabricate false precision
(it would make the ambiguity set look smaller than the evidence supports); that
is explicitly forbidden.

Every use of incomplete reference is surfaced in `ObservationResult.loss_record`
and in dedicated metadata fields (§6.2), so downstream analyses can exclude or
weight these observations.

---

## 5. Resolution modes expressed on the engine

All M2 sequence modes are *configurations* of the same engine (spec §3, §4),
differing only in `interrogated_features` / masks and flags — not separate code
paths:

| Mode | Interrogated region | Notes |
|------|--------------------|-------|
| `second_field` (sequence) | assay-defined exons | Replaces idealized truncation; may yield `/` or `\|` ambiguity (§7). |
| `intermediate` | discriminatory probes / partial exons / explicit signature set | Ambiguity is whatever the probes cannot split; **not** a MAC lookup (spec §4/§5). |
| `exon2_3_classI` | exon 2 + exon 3 | Class I preset. |
| `exon2_classII` | exon 2 | Class II preset. |
| `full_gene` | all features with `PRESENT` sequence | Resolves to the limit of complete reference; need **not** resolve all nomenclature distinctions (spec §5) — e.g. intronic-only differences with absent intron data remain ambiguous. |

The idealized M1 modes (`allele`, `first_field`, `second_field`-idealized)
remain available and unchanged; `assay_class = idealized` still routes to the M1
per-copy path. `assay_class ∈ {sso, sbt}` routes to the M2 engine. This keeps
the M1 guarantee that idealized transforms are never presented as empirical
(spec §12).

---

## 6. GL String output and representation fidelity

### 6.1 Normal output

The engine's per-locus `LocusGroup` flows into the existing serializer. Because
genotype ambiguity is modeled as `LocusGroup.alternatives` (the `|` layer that
already exists in the M1 AST), **no GL String grammar extension is required**
for the common cases. Allele ambiguity uses the existing `AlleleList` (`/`).

### 6.2 Metadata additions (backward-compatible)

`ObservationMetadata` / `LocusObservation` gain sequence-specific fields. These
are additive; M1 results are unaffected.

```
LocusObservation (extended):
    interrogated_features: tuple[FeatureKey, ...]
    candidate_genotype_count: int
    reduction: "unambiguous" | "allele_factored" | "genotype_ambiguous"
               | "superset_fallback"
    uses_incomplete_reference: bool
    indeterminate_alleles: tuple[str, ...]

ObservationResult (extended):
    representation_limitations: tuple[RepresentationLimitation, ...]
```

`error_category` stays `NONE` for all of M2: assay ambiguity and information
loss are **not** mistyping (spec §7). M3 is the only layer that sets non-`NONE`.

### 6.3 When the representation cannot be faithful

Two known cases where a naive serialization would lie:

- **Cross-locus constraint from a shared amplicon** (rare; e.g. an assay whose
  signature couples two loci). GL String cannot encode "these two `|` sets are
  jointly constrained." The engine must either keep the loci in one `LocusGroup`
  with explicit joint alternatives (if they share a block) or emit a
  `RepresentationLimitation` and fall back to independent per-locus supersets.
- **Partial intralocus phase** that is neither full cis/trans resolution nor
  full ambiguity. If the grammar in use cannot express the partial constraint,
  emit the safe superset and flag it.

Rule (spec §8): **never silently flatten or narrow.** A superset fallback is
acceptable and flagged; a narrowing is never acceptable.

### 6.4 Phase regrouping must not invent cross-locus associations

Correction (review item 6). Per-locus candidate enumeration (§3.2) is done
**independently** at each locus. The regrouping step that reassembles loci into
a `Diplotype` must respect exactly the phase the assay observed, with two rules:

**Rule A — no joint phase observed ⇒ independent groups, no cross-product.**
If the assay does not jointly phase loci B and C, the engine emits the two
per-locus candidate sets as **separate** `^`-delimited groups. It must **not**
form the cross-product of B-haplotypes with C-haplotypes, because that would
assert an interlocus association the assay never saw.

**Rule B — joint phase observed ⇒ evaluate candidates against joint phase.**
If the assay jointly phases B~C (e.g. long-range SBT over a B–C amplicon, or a
declared phase group backed by a phase-resolving assay), candidate *haplotype
pairings* must be enumerated and compared against the **jointly observed phase**,
not reconstructed by re-attaching phase to independently chosen per-locus pairs.

#### Synthetic example (verified)

Truth, phased: haplotypes `B*A~C*A` and `B*B~C*B` (cis: B\*A with C\*A, B\*B with
C\*B). Suppose each locus is individually unambiguous: B observation ⇒ exactly
`{B*A, B*B}`, C observation ⇒ exactly `{C*A, C*B}`.

- **Assay does NOT phase B~C** (`resolves_phase` false, no joint observation):
  ```
  B*A+B*B ^ C*A+C*B
  ```
  This asserts no interlocus relationship — correct. Note it is *consistent with
  both* cis assignments `{B*A~C*A, B*B~C*B}` (truth) and `{B*A~C*B, B*B~C*A}`.
  That is the honest state of knowledge: the assay could not tell them apart.
  **The engine must stop here and not re-attach phase.** Re-attaching the truth
  phase at this point would invent an association (`B*A~C*A`) the assay did not
  observe — the exact error this rule forbids.

- **Assay DOES phase B~C** (joint observation of the two haplotypes):
  the two cis assignments are **distinct genotypes**
  ```
  truth:  B*A~C*A + B*B~C*B
  alt:    B*A~C*B + B*B~C*A
  ```
  Joint evaluation compares each candidate *haplotype pairing* against the
  observed phase. Only the pairing(s) the assay actually resolved are kept. If
  the assay resolved truth's cis phase, the output is the single phased block
  ```
  B*A~C*A + B*B~C*B
  ```
  and the alt pairing is excluded — because it was observed to be incompatible,
  not merely unobserved. If the assay phased the region but could not break the
  cis/trans symmetry (e.g. both haplotypes share the phased columns), **both**
  pairings are emitted as genotype ambiguity within one phased group:
  ```
  B*A~C*A+B*B~C*B | B*A~C*B+B*B~C*A
  ```
  This `|` lives *inside* a single `LocusGroup` spanning {B,C} (it is genotype
  ambiguity over phased haplotypes), which the M1 AST already expresses.

The distinction between the unphased result (`B*A+B*B ^ C*A+C*B`) and the phased
ambiguous result (`…~… | …~…`) is the whole point: the first asserts no linkage;
the second asserts linkage exists but is two-way ambiguous. The engine must
produce whichever the assay's phase observation justifies, and never upgrade the
first into the second by re-attaching known truth phase.

**Implementation note.** This means joint-phase loci cannot be enumerated purely
per-locus-then-regroup. When a declared phase group is backed by a phase-
resolving assay, the engine enumerates candidate *haplotype pairings* over the
joint interrogated columns of the grouped loci (the §3.2 algorithm applied to
the concatenated multi-locus signature), and only then reduces to GL form.
Loci outside any phase-resolving group retain the independent per-locus path.

#### 6.4.1 Worked example — per-locus ambiguities correct, cross-product wrong

This is the case review item 6 asks for: **each locus is individually ambiguous
and individually correct, yet the unrestricted combination of the two per-locus
ambiguities introduces a multilocus diplotype that the joint phase observation
excludes.** (Verified numerically.)

Setup. A B~C phase-resolving assay (long-range SBT across a B–C molecule). Truth
is phased into two molecules:

```
molecule 1 (hap1):  B*A  ~ C*1
molecule 2 (hap2):  B*2  ~ C*B
```

The assay cannot resolve certain intra-signature distinctions, so **each locus
is genuinely ambiguous**:

- at **B**: `B*A` is indistinguishable from `B*Aq` (same interrogated B
  columns). The B genotype is `B*A/B*Aq + B*2`.
- at **C**: `C*B` is indistinguishable from `C*Bq`. The C genotype is
  `C*1 + C*B/C*Bq`.

Both per-locus genotypes are **correct** statements of what each locus assay
alone resolves.

**Unrestricted combination (WRONG).** Treating the two per-locus genotypes as
independent and re-attaching phase arbitrarily yields a B-copy paired with
either C-copy. That fabricates, among others:

```
B*A ~ C*B  +  B*2 ~ C*1          <-- NOT observed on any molecule
```

Here `C*B` is placed in cis with `B*A`, but the assay read `C*B/C*Bq` **on the
same molecule as `B*2`**, and `C*1` **on the same molecule as `B*A/B*Aq`**. The
fabricated diplotype asserts a haplotype (`B*A~C*B`) that no molecule carried —
an invented cross-locus allele association.

**Joint-phase evaluation (CORRECT).** Because the assay resolved B~C phase, the
engine enumerates candidate *haplotype pairings* against the two observed
molecules. Molecule 1's signature admits `{B*A~C*1, B*Aq~C*1}`; molecule 2's
admits `{B*2~C*B, B*2~C*Bq}`. The candidate diplotype set is the product **of
molecules, not of loci** — exactly 4 phased diplotypes, each a valid
`LocusGroup` over {B,C}:

```
B*A~C*1  + B*2~C*B    |  B*A~C*1  + B*2~C*Bq   |
B*Aq~C*1 + B*2~C*B    |  B*Aq~C*1 + B*2~C*Bq
```

The fabricated `B*A~C*B + B*2~C*1` is **absent** (verified: it is not in the
joint candidate set). The contrast in counts is the teaching point:

| Enumeration | Diplotypes | Correct? |
|-------------|-----------|----------|
| product of per-locus genotypes, phase re-attached arbitrarily | 8 (4 cis + 4 fabricated trans) | **No** — invents cross-locus haplotypes |
| product of observed molecules (joint phase) | 4 | **Yes** |

If instead the assay did **not** resolve B~C phase, the honest output is the two
independent groups `B*A/B*Aq+B*2 ^ C*1+C*B/C*Bq`, which asserts *no* interlocus
pairing at all — neither the 4 cis nor the 4 trans — and must not be "upgraded"
by re-attaching truth phase (§6.4 Rule A).

---

## 7. Worked examples: why independent allele matching is wrong

These are the crux of the milestone. In each, "independent matching" means:
find the candidate set for `A1` alone, find the candidate set for `A2` alone,
then take their cross product and serialize as `(/-list) + (/-list)`. The engine
instead enumerates **pairs** that reproduce the combined observation (§3.2) and
only factorizes to `/` when the candidate set is exactly a Cartesian product
(§3.3 step 2).

Notation below uses toy two-column signatures `[c1 c2]` for clarity; real
columns are alignment positions over interrogated exons.

### 7.1 Independent matching invents a non-existent homozygote

Truth at a locus: heterozygous `A*01 + A*02`, with interrogated signatures

```
A*01 = [G A]      A*02 = [T C]
```

Suppose two other alleles exist with

```
A*03 = [G C]      A*04 = [T A]
```

Combined observation (unphased, set semantics):

```
O = [ {G,T}  {A,C} ]
```

**Independent matching** would ask "what matches A*01's `[G A]`?" → {A*01}, and
"what matches A*02's `[T C]`?" → {A*02}, giving the single genotype
`A*01 + A*02`. That looks fine here — but now suppose the assay only reads
**column 1** (`interrogated = {c1}`):

```
signatures over c1:  A*01=G  A*03=G   |   A*02=T  A*04=T
O = { {G, T} }
```

Independent matching: `A*01` matches {A*01, A*03}; `A*02` matches {A*02, A*04};
cross product → `A*01/A*03 + A*02/A*04`. This is actually the correct Cartesian
set here, so `/` factorization is valid. The point of 7.1 is the *contrast* with
7.2, where the Cartesian assumption breaks.

### 7.2 The real failure: cross-pairing the `/` product asserts excluded genotypes

Keep the single-column assay, but change the allele inventory so the Cartesian
product contains pairs the observation **excludes**. Interrogated columns = both
`c1` and `c2`, assay is **not** phase-resolving, set semantics:

```
A*01 = [G A]   (true)
A*02 = [T C]   (true)
A*05 = [G C]
A*06 = [T A]
```

Combined true observation:

```
O = [ {G,T}  {A,C} ]
```

Which **pairs** reproduce `O`? Enumerate (§3.2):

- `{A*01, A*02}` = `[{G,T} {A,C}]` ✓
- `{A*05, A*06}` = `[{G,T} {C,A}]` ✓   (G&C from A*05, T&A from A*06)
- `{A*01, A*06}` = `[{G,T} {A,A}]` = `[{G,T} {A}]` ✗ (c2 not heterozygous)
- `{A*02, A*05}` = `[{T,G} {C,C}]` ✗
- homozygotes ✗

So the **correct** genotype ambiguity is:

```
A*01+A*02 | A*05+A*06
```

**Independent matching** would compute: alleles matching the true A*01 at both
columns = {A*01}; matching A*02 = {A*02}. That yields only `A*01+A*02` and
**misses** the genuinely indistinguishable `A*05+A*06` — an under-ambiguous,
falsely precise result.

Worse, a different naive variant — "any allele consistent with *some* observed
base at each column" — gives per-chromosome candidate set {A*01, A*02, A*05,
A*06} on both sides, and the `/`-cross-product

```
A*01/A*02/A*05/A*06 + A*01/A*02/A*05/A*06
```

which **asserts** excluded genotypes like `A*01+A*06` (shown ✗ above). This is
the canonical independent-matching error: using `/` (allele ambiguity, which
distributes over the pairing) to represent what is really a constrained **set of
pairs** (`|`, genotype ambiguity). The two operators are not interchangeable,
and conflating them is precisely the bug §4 of the spec warns against.

The engine avoids both failure modes because step §3.3-2 only emits `/` when the
candidate set equals the full Cartesian product; here it does not
(`A*01+A*06` is absent), so it emits the `|` form.

### 7.3 Set vs multiset: the only regime where they differ is out of M2 scope

This example states its assumptions explicitly and shows that for a **clean
diploid** a set and a multiset carry identical information — so the removed
`zygosity_aware` flag had no effect on a valid M2 observation — and isolates the
*only* regime where a multiset would differ (a >2-chromosome mixture), which is
M3.

**Assumptions:** exactly two chromosomes contribute (clean diploid); the assay
reads both columns completely (no dropout); reference is PRESENT for all alleles
shown.

Truth: homozygous `A*01 + A*01`, with `A*01 = [G A]`. A second allele
`A*07 = [G A]` is identical over the interrogated columns (differs only in an
intron the assay does not read).

- **Observed column sets (the M2 model):**
  `O = [ {G}  {A} ]`. Each column is a singleton, which — because the
  observation is complete (§2.3.1) — is a *confirmed homozygous* reading.
- **As a multiset it would be** `[ {G:2} {A:2} ]`. For a diploid this carries
  **no extra information**: a singleton set already implies both chromosomes
  equal that base (there is no third copy whose count could vary). So set and
  multiset agree.
- **Candidate pairs** (identical interrogated sequence to truth):
  `{A*01,A*01}`, `{A*01,A*07}`, `{A*07,A*07}`. This set equals the Cartesian
  product `{A*01,A*07} × {A*01,A*07}`, so the exact reduction (§3.3) is
  `A*01/A*07 + A*01/A*07`. (Verified: 3 unordered pairs = the product.)

The candidate set above is **independent of any set/multiset choice** — the
reviewer's point. The previous draft's claim that "set semantics cannot tell
homozygous from heterozygous" was the error; for a complete diploid observation
the set `{G}` *is* the homozygous call.

**Where a multiset would genuinely differ — and why it is M3.** Suppose instead
that **three** chromosomes contribute (sample mixture / contamination), with
bases `G, G, T` at a column:

- set: `{G, T}` (looks like an ordinary heterozygote);
- multiset / dosage: `{G:2, T:1}` (reveals the 2:1 imbalance → evidence of a
  mixture).

Only here does dosage add information, and a 3-contributor mixture is precisely
the contamination/mixture mechanism deferred to **M3** (spec §7). M2 models the
clean, complete diploid, where `resolves_dosage` is inert (§2.3.1). This
example is the analytic justification for removing the M2 flag.

### 7.4 Missing reference keeps truth in the set without faking identity

Truth `A*01 + A*08`. `A*08` has `ABSENT` exon-3 data; assay interrogates exon 2
+ exon 3. Over exon 2 `A*01` and `A*08` differ (INCOMPATIBLE) → they are
distinguishable there, truth pair stands. But a third allele `A*09` matches
`A*01` on exon 2 and is `INDETERMINATE` vs `A*08` on exon 3 (A*08 absent). The
pair `{A*01, A*09}` is a candidate **only because exon-3 is indeterminate**, so
it is included, the result is flagged `uses_incomplete_reference = True`, and
`A*08, A*09` are listed in `indeterminate_alleles`. The engine never records the
indeterminate exon-3 columns as a match (§4 rule 2).

---

## 8. Algorithm summary and performance

Per locus:

1. Build `U` = alleles with usable sequence over interrogated columns (§3.1, §4).
2. Compute `signature(u)` for `u ∈ U`; bucket alleles by identical signature.
3. Compute true `O = combined_signature(A1, A2)` (§2.3).
4. Candidate pairs = pairs of buckets whose combined signature equals `O`
   (three-valued comparison, §4). Bucketing makes this `O(k^2)` in the number of
   **distinct** signatures `k`, not `O(|U|^2)`.
5. Reduce to the tightest faithful GL form (§3.3); flag any superset fallback.
6. Emit `LocusGroup` + extended metadata (§6.2).

Determinism: canonical sorting of alleles, pairs, and alternatives guarantees
stable, round-trippable output (QA §10.9) and reproducibility with a fixed seed
(no randomness in M2; QA §10.8 holds trivially).

---

## 9. Serology (structural sketch, for completeness)

Serology is **not** a truncated DNA name (spec §5). It is a versioned
`allele -> {serologic specificities}` mapping with broad/split relationships and
explicit "unresolved"/"blank" assignments. The engine reuses the same
pair-enumeration skeleton: the "signature" of an allele is its serologic
specificity set; the combined observation is the union of specificities across
the two alleles with broad/split collapsing applied by the versioned mapping;
candidate pairs are those reproducing the observed specificity multiset. Output
uses serologic names (e.g. `HLA-B7`, `HLA-Bw4`) in the allele positions.

The full serology design — reference artifact, broad/split lattice, blank vs
unknown handling, worked examples, and QA — is specified in
`docs/design/m2-serology.md`. It plugs into the shared genotype-ambiguity engine
as a second `ObservationBackend` (decision D6) with no change to §3.2–§3.3.

---

## 10. QA mapping (what M2 tests must assert)

| Spec §10 item | M2 obligation |
|---|---|
| 1 error-free obs compatible with truth | truth pair ∈ candidate set, always |
| 4 unphased discards B~C | unchanged from M1; verify engine adds no phase (§6.4 Rule A) |
| 4b multilocus phase not invented | §6.4.1: joint-phase candidates = product of observed molecules, not of loci; no fabricated cross-locus haplotype |
| 5 ambiguity set contains truth | §3.2 guarantees `{A1,A2} ∈ candidates` |
| 6 errors can be incompatible | out of scope (M3); M2 never incompatible |
| 8 reproducible with seed | deterministic; verify byte-identical output |
| 9 valid round-trippable GL | `/` vs `\|` reduction preserves round-trip |
| 10 reference versions recorded | `database_version` threaded into metadata |
| challenge: incomplete reference | §4; verify no indeterminate→match promotion |
| challenge: null alleles | null expression carried as sequence feature; a null from a coding change in the interrogated region is distinguishable, one from outside is not — verify both |
| challenge: cis/trans genotype ambiguity | §7.2 case produces `\|`, never cross-`/` |
| challenge: DRB3/4/5 | locus-presence (gene may be absent) modeled as "no genotype at that locus", distinct from "typed, homozygous" |

Three regression tests are mandatory and encode the thesis of this milestone:

- **`test_no_independent_cross_product`** — the §7.2 inventory must yield
  `A*01+A*02 | A*05+A*06` and must **not** yield any `/`-cross-product
  containing `A*01+A*06`.
- **`test_truth_always_in_candidates`** — randomized over the reference set.
- **`test_incomplete_reference_not_treated_as_match`** — the §7.4 case.
- **`test_unphased_regroup_asserts_no_interlocus_pairing`** (§6.4 Rule A) — for a
  non-phase-resolving B+C strategy over phased truth, the result is two
  independent `^` groups and contains **no** `~` linking B to C; assert the
  engine never re-attaches truth phase.
- **`test_joint_phase_product_of_molecules_not_loci`** (§6.4.1) — the §6.4.1
  inventory under a B~C phase-resolving assay yields exactly the **4** phased
  diplotypes (product of observed molecules) and must **not** contain the
  fabricated `B*A~C*B + B*2~C*1` (nor any of the 4 trans recombinants).
- **`test_phase_regroup_never_invents_haplotype`** — property test: for random
  phased truth and a non-phase-resolving strategy, every haplotype appearing in
  the output must be a sub-haplotype actually present in truth (no invented
  cross-locus cis association).

---

## 11. Open questions for review — RESOLVED

All five questions below have been resolved in
`docs/design/m2-decisions.md` (decisions D1–D5), with two further decisions
(D6 backend seam, D7 null semantics) recorded there. Serology is specified in
`docs/design/m2-serology.md`. Summary of resolutions:

1. **Alignment source of truth** → precomputed, checksummed per-locus column
   model artifact (D1).
2. **Column model for indels** → `gap` and `indeterminate` are distinct symbols;
   sequenced deletions are evidence, missing data is not (D2).
3. **Multiset vs set default** → `zygosity_aware` is required per sequence assay,
   preset per assay physics (SSO/SBT = set, NGS/long-read = multiset) (D3).
4. **Superset fallback policy** → `representation_policy ∈ {superset(default),
   strict, annotate_only}`, never narrowing (D4).
5. **Serology companion** → serology stays in M2 as a second `ObservationBackend`
   (D5/D6), specified in `docs/design/m2-serology.md`.

Nothing above or in the companions is implemented; these documents are the
agreed basis for the M2 implementation once review signs off.
