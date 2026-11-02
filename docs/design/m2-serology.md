# M2 Serology Companion Design

**Status:** design only. Nothing here is implemented.
**Companion to:** `docs/design/m2-sequence-ambiguity.md` (main M2 design) and
`docs/design/m2-decisions.md` (decisions D1–D7). This document specifies the
**serology observation backend** referenced by D5/D6/D7.
**Scope:** spec §5 (Serology), §4 (realistic ambiguity), §8 (GL String), §10 QA.

> **Principle (spec §5):** a serologic result is **not** a truncated DNA allele
> name. It is the antigen the antiserum recognizes, via a versioned
> allele→specificity mapping with broad/split structure and explicit blanks.

The serology backend reuses the shared genotype-ambiguity algorithm from the
main design (§3.2–§3.3) unchanged. It only supplies a different **signature**
(serologic specificities instead of base columns) and a different **reference
artifact**. This document specifies that signature, the broad/split lattice, the
reference data, blank/unresolved handling, worked examples, GL String output,
and the QA obligations.

---

## 1. Why serology is a backend, not a name transform

A naive serology model truncates `HLA-B*07:02` → `B7`. That is wrong in both
directions:

- Many alleles of one DNA family map to **different** serologic specificities
  (and vice versa); the relationship is a curated many-to-many mapping, not a
  prefix.
- Some alleles have **no** serologic equivalent (blank), or only a **broad**
  assignment with the split unresolved.
- A serologic *typing* observes antigens on cells from **both** chromosomes at
  once — so, exactly as in the sequence backend (main §2.3), the observation is
  a property of the **pair**, and genotype ambiguity must be computed over
  pairs, not by mapping each allele independently.

So serology plugs into the same `ObservationBackend` seam (decision D6):

```
SerologyBackend.signature(allele, config)          -> SpecificitySet
SerologyBackend.combined_observation(a1, a2, config) -> SpecificityObservation
```

The shared engine enumerates DNA-allele pairs whose `combined_observation`
reproduces truth's and reduces to the **compatible-DNA GL String** (§7) — exact
`/`-vs-`|` discipline — alongside a separate structured **antigen report**.
Identical control flow to the sequence backend; only the signature and output
packaging differ.

---

## 2. Reference data: the versioned serology artifact

Parallel to the column-model artifact (decision D1), serology uses a pinned,
checksummed artifact per database version. Source of truth is the
IPD-IMGT/HLA serology relationships (the `rel_dna_ser.txt` assignment table plus
the WHO broad/split relationships).

```
SerologyModel (per database_version):
    database_version: str                 # e.g. "IPD-IMGT/HLA 3.55.0"
    assignment: Mapping[Allele, Assignment]
    lattice: BroadSplitLattice
    checksum: str

Assignment:
    specificities: frozenset[Specificity]   # may be empty (blank)
    confidence: ASSIGNED | POSSIBLE | UNKNOWN | NOT_EXPRESSED
    #   ASSIGNED  '+'  in rel_dna_ser   (expressed, confirmed)
    #   POSSIBLE  '?'                    (expected but unconfirmed)
    #   UNKNOWN   blank cell             (no assignment on record)
    #   NOT_EXPRESSED  '0'               (null / not expressed -> blank antigen)
```

`confidence` is retained because it changes candidate inclusion (§5): a
`POSSIBLE` specificity must not be treated as a hard observation, and an
`UNKNOWN` must never be treated as "matches anything."

The `database_version` + `checksum` are threaded into
`ObservationResult.metadata.reference_databases` (QA §10.10), exactly as for
sequence.

---

## 3. The broad/split lattice

Serologic specificities form a shallow hierarchy: a **broad** antigen (e.g.
`B15`) splits into **splits** (e.g. `B62`, `B63`, `B75`, `B76`, …). Some
antisera resolve the split; some report only the broad.

```
BroadSplitLattice:
    broad_of: Mapping[Specificity, Specificity]   # split -> its broad (or self)
    splits_of: Mapping[Specificity, frozenset]    # broad -> its splits
    associated: Mapping[Specificity, frozenset]   # e.g. Bw4/Bw6 public epitopes
```

The lattice is applied by a **resolution level** declared on the assay config:

- `split`  — the assay resolves splits; observation uses split-level
  specificities.
- `broad`  — the assay reports only broads; every specificity is projected up to
  `broad_of` before comparison. This is where historical serology loses
  information and genotype ambiguity appears (worked in §6.2).

`associated` public epitopes (Bw4/Bw6, and the C locus `Bw` equivalents) are
modeled as **additional** specificities an allele carries, not as a separate
axis — so an antiserum panel that reports Bw4 contributes a `Bw4` token to the
signature just like any other specificity.

---

## 4. Signature and combined observation

**Signature of one allele** under the serology backend at resolution `R`:

```
signature(allele, R) =
    project( assignment[allele].specificities , R )     # split or broad level
    filtered by confidence policy (§5)
```

**Combined observation** of the true pair `{A1, A2}` is the **union** of the two
alleles' projected specificity sets (set semantics by default, because a
serologic panel reports antigen presence, not dosage):

```
combined_observation(A1, A2, R) =
    signature(A1, R)  ∪  signature(A2, R)
```

The serologic observation is a **set** of antigen tokens (presence/absence): a
serologic panel reports which antigens are present, not their dosage. Unlike the
sequence backend — where for a clean diploid the base *set* already encodes
zygosity (main §2.3.1) — a serologic antigen set genuinely **cannot** encode
zygosity, because one antigen can be contributed by one or both chromosomes and
the panel reports only "present." There is no multiset/dosage variant for
standard serology; the removed `zygosity_aware` flag (decision D3) does not
apply here.

> Consequence: serology **cannot distinguish homozygous from heterozygous** for
> the same antigen. `A1,A2` both `B7` → observation `{B7}`, indistinguishable
> from a single-antigen case. This is honest ambiguity that feeds the candidate
> enumeration (§6.3), **not** a false-homozygosity *error* (that is M3). It is a
> genuine property of the serologic observation, not an assay-model choice.

---

## 5. Confidence and blank/unresolved handling

The confidence field gates how a specificity participates, preventing the two
classic faking errors (never invent an antigen, never treat "no data" as
evidence):

| Confidence | In truth's observation | In a candidate's signature |
|------------|------------------------|-----------------------------|
| `ASSIGNED` | contributes the antigen | must be able to produce it |
| `POSSIBLE` | contributes, flagged `uses_possible_assignment` | may produce it (not required) |
| `UNKNOWN` | contributes **nothing**; allele's antigen is indeterminate | the allele is INDETERMINATE for matching — included only as a flagged candidate, never as a confirmed match |
| `NOT_EXPRESSED` (null) | contributes **blank** (no antigen) | must likewise be blank |

**Blank handling.** A blank specificity set is a real observation outcome
(serologically "not detected"), **distinct** from UNKNOWN (no data on record).
The engine represents:

- blank (`NOT_EXPRESSED` or genuinely no antigen) → the allele adds nothing to
  the union, and candidate pairs must reproduce that absence;
- UNKNOWN → the allele is indeterminate; candidates relying on it are flagged
  (`uses_unknown_assignment = True`) and recorded, exactly parallel to
  incomplete-reference handling for sequence (main §4, decision D7).

This is the serologic form of "handle missing reference explicitly; do not
silently treat unavailable data as identity" (spec §4).

---

## 6. Worked examples

Verified by hand-computation; the arithmetic matches the engine's pair
enumeration (main §3.2). Toy specificity map:

```
B*07:02 -> {B7}       B*08:01 -> {B8}
B*15:01 -> {B62}      B*15:02 -> {B75}     B*15:03 -> {B72}
B*15:99 -> {B15}      (broad only; split unresolved on record)
lattice: B62,B63,B75,B72,B15 all broad_of -> B15
```

### 6.1 Split-level assay: unambiguous

Truth `B*15:01 + B*07:02`, resolution `split`. The result has **two separate
parts** (review item 3):

```
Antigen report (structured, NOT a GL String):
    specificities observed = { B62, B7 }        # split-level antigens
    resolution_level = split
    public_epitopes = { Bw6 }                    # if the panel typed Bw4/Bw6 (see §6.5)

Compatible DNA genotype set (a valid GL String of DNA alleles):
    HLA-B*07:02+HLA-B*15:01
```

Only one DNA pair yields `{B62, B7}` here, so the compatible-genotype GL String
is unambiguous. **`B62` and `B7` never appear in the GL String** — they are not
nomenclature-level allele identifiers; they live only in the antigen report.

### 6.2 Broad-only (historical) assay: genotype ambiguity appears

Same truth, resolution `broad`. Every B15-group split projects to `B15`:

```
O = {B15} ∪ {B7} = {B15, B7}     # broad-level antigens only
candidate DNA pairs reproducing O =
    { {B*07:02, B*15:01}, {B*07:02, B*15:02},
      {B*07:02, B*15:03}, {B*07:02, B*15:99} }
```

Two separate parts again:

```
Antigen report (structured):
    specificities observed = { B15, B7 }         # BROAD only — never B62/B75/...
    resolution_level = broad

Compatible DNA genotype set (GL String, DNA alleles only):
    HLA-B*07:02+HLA-B*15:01/HLA-B*15:02/HLA-B*15:03/HLA-B*15:99
```

Key corrections (review item 3):

- **Never report split-level information from a broad-only assay.** The antigen
  report lists `{B15, B7}` — broads. It does **not** list `B62/B75/B72`; the
  assay did not resolve the split, so claiming a split would fabricate an
  observation. The split identities survive only implicitly, as the *DNA
  alleles* in the compatible set, which is a statement about "alleles consistent
  with the broad antigen," not "splits observed."
- The compatible-genotype GL String contains **only DNA nomenclature alleles**
  joined per the exact `/`-vs-`|` discipline (main §3.3). Here the set factors
  exactly as a Cartesian product (left fixed `B*07:02`, right any B15-group
  allele), so `/` is valid; a non-factorable set would use `|` (see §6.6).
- The engine emits the antigen report and the compatible-DNA GL String as
  distinct fields; it never serializes antigen tokens as GL String alleles.

### 6.3 Homozygous-looking: serology cannot see the second antigen

Truth `B*07:02 + B*08:01` vs truth `B*07:02 + B*07:02`:

```
het:  O = {B7, B8}
hom:  O = {B7}
```

Now consider truth `B*07:02 + B*07:02` (`O = {B7}`). Candidate pairs reproducing
`{B7}` include **every** pair whose union is exactly `{B7}`: `{B*07:02,
B*07:02}` and any other allele also mapping solely to `B7` (homozygous or
heterozygous among them). Serology genuinely cannot assert homozygosity; the
engine emits the honest genotype-ambiguity set and flags it. This is **not** the
M3 false-homozygosity error — truth remains in the set (QA §10.5).

### 6.4 Null allele → blank (decision D7)

Truth `B*07:02 + B*15:01:01:02N` (a null). The null is `NOT_EXPRESSED` → blank.

```
O = {B7} ∪ {} = {B7}
```

The observation is serologically identical to a `B7` homozygote and to any
`B7 + <null>` pair. Candidate set therefore includes `{B*07:02, B*07:02}` and
`{B*07:02, <any null or blank allele>}`. This is exactly how serology's
blindness to a null is modeled honestly (D7): the engine does **not** infer the
null from the `N` suffix; it reads `NOT_EXPRESSED` from the serology artifact.

### 6.5 Public epitopes Bw4/Bw6 are a separate axis of the antigen report

Bw4 and Bw6 are **public epitopes** carried by many B (and some A, C) alleles,
reported by a panel independently of the broad/split specificity. They are
**not** broad/split antigens and **not** DNA alleles. They occupy their own
field in the antigen report:

```
Truth B*07:02 + B*08:01   (B*07:02 -> Bw6, B*08:01 -> Bw6)
Antigen report:
    specificities = { B7, B8 }
    public_epitopes = { Bw6 }            # both carry Bw6; homozygous-looking at this axis
```

A pair carrying one Bw4 and one Bw6 allele reports `{Bw4, Bw6}`. The epitope
axis participates in candidate comparison exactly like a specificity (a pair is
incompatible if it cannot reproduce the observed epitope set), but it is
recorded and reported separately so a consumer can see "Bw4/Bw6 typed" distinct
from "B-antigen typed." Epitope tokens never enter the compatible-DNA GL String.

### 6.6 The naive union cross-product asserts excluded genotypes

The serology analogue of main §7.2: do **not** build the compatible-DNA set by
taking every allele consistent with *some* observed antigen and crossing it with
itself. Toy artifact (broad level) with one **cross-reactive** allele:

```
a = B*07:02 -> {B7}
b = B*08:01 -> {B15}
c = B*15:16 -> {B7, B15}      # broad cross-reactive (carries both)
```

Observed broad antigens `O = {B7, B15}`. Enumerating unordered DNA pairs whose
antigen union equals `O` (verified):

```
candidate pairs = { a+b, a+c, b+c, c+c }
```

A **naive union cross-product** — "alleles touching B7 or B15" = `{a,b,c}`,
crossed with itself — would emit `a/b/c + a/b/c`, which asserts `a+a` (`{B7}`
only) and `b+b` (`{B15}` only). **Both are excluded** by the observation
(verified: `naive − candidates = {a+a, b+b}`). That is the independent-matching
error.

The exact candidate set *does* factor here as `{a,c} × {b,c}`, so the exact,
lossless GL form is:

```
HLA-B*07:02/HLA-B*15:16 + HLA-B*08:01/HLA-B*15:16
```

The engine must emit this exact factorization (or, when a set does not factor,
exact `|` enumeration — main §3.3, default `exact`), **never** the naive union
cross-product. The rule is identical to the sequence backend: `/` only when it
equals the true candidate set; otherwise exact `|`.

---

## 7. Output contract: antigen report vs compatible DNA set

Correction (review item 3). A serologic observation produces **two distinct
outputs**, never conflated:

**(1) The antigen report** — a structured object, *not* a GL String:

```
SerologicAntigenReport (per locus):
    specificities: frozenset[AntigenToken]       # broad or split, per resolution_level
    resolution_level: BROAD | SPLIT
    public_epitopes: frozenset[EpitopeToken]     # Bw4/Bw6 etc., separate axis
    confidence_flags: {uses_possible_assignment, uses_unknown_assignment, blank_present}
```

`AntigenToken` (`B7`, `B62`, `B15`) and `EpitopeToken` (`Bw4`, `Bw6`) are
**serologic identifiers, not HLA nomenclature alleles.** They are deliberately a
different type from `Allele` and **must never be placed in a GL String**. A
broad-only assay sets `resolution_level = BROAD` and lists only broad tokens —
split tokens are never emitted (§6.2).

**(2) The compatible DNA genotype set** — a valid GL String over DNA
nomenclature alleles:

- It is the set of DNA genotypes consistent with the antigen report, serialized
  with the exact `/`-vs-`|` discipline (main §3.3): `/` only when the set is an
  exact Cartesian product (§6.2), otherwise exact `|` enumeration (§6.6).
- It contains only `HLA-<locus>*<fields>` tokens — the same `Allele` type the M1
  parser accepts. No antigen or epitope token ever appears here.
- Representation policy (decision D4) applies: **exact by default**; a superset
  only on explicit `representation_policy = superset`, flagged in metadata.

**No GL String grammar change is needed** — but for the correct reason: the GL
String carries only DNA alleles (which the M1 grammar already handles), and the
antigen report is carried *outside* the GL String as structured metadata. The
earlier claim that `HLA-B7`/`HLA-Bw4` are "valid GL String allele tokens" is
withdrawn; they are not nomenclature alleles and are not serialized into GL.

`report_vocabulary` is removed: the result always carries both parts (the
antigen report and the compatible-DNA GL String), so a consumer chooses which to
read rather than the engine choosing one vocabulary to emit.

Representation limitations (cross-constraint, inexpressible partial ambiguity)
follow decision D4 policy (`exact` default / `superset` / `annotate_only`).

---

## 8. QA obligations (serology-specific, extending main §10)

| Item | Serology obligation |
|------|---------------------|
| truth in candidate set (§10.5) | union-reproduction guarantees `{A1,A2}` ∈ candidates |
| broad/split correctness | §6.1 unambiguous at split; §6.2 four-way ambiguous at broad |
| blank ≠ unknown | §5 table; a `NOT_EXPRESSED` null yields blank (§6.4), an `UNKNOWN` allele is flagged indeterminate, never a confirmed match |
| homozygous vs heterozygous | §6.3 — serology cannot assert homozygosity; truth stays in set; not labeled error |
| not a name truncation | regression test: an allele whose serologic specificity is **not** a prefix of its DNA name (e.g. `B*15:01 -> B62`) must map via the artifact, never by truncation |
| reference version recorded (§10.10) | serology artifact `database_version` + `checksum` in metadata |
| cis/trans genotype ambiguity (§10) | §6.2 produces the correct `/`-factorable case; a constructed non-factorable case must stay `|` (no cross-`/`), same thesis as sequence §7.2 |
| DRB3/4/5 locus presence | a locus with no serologic equivalent on an absent gene is "no genotype," distinct from blank-but-present |

Mandatory regression tests mirroring the sequence thesis:

- **`test_serology_broad_split_ambiguity`** — §6.2: antigen report is `{B15,B7}`
  (broad only, no split tokens), and the compatible-DNA GL String is
  `HLA-B*07:02+HLA-B*15:01/HLA-B*15:02/HLA-B*15:03/HLA-B*15:99`.
- **`test_serology_tokens_never_in_gl_string`** — assert no antigen/epitope token
  (`B7`, `B62`, `Bw4`) ever appears in the compatible-DNA GL String; they live
  only in the structured antigen report.
- **`test_serology_broad_only_never_emits_split`** — a broad-only assay's antigen
  report must contain no split-level token under any inventory.
- **`test_serology_is_not_name_truncation`** — `B*15:01` serologic specificity is
  `B62` (via artifact), never `B15` from the name.
- **`test_serology_null_is_blank_not_inferred`** — §6.4; the null's blank comes
  from `NOT_EXPRESSED`, and removing the artifact entry makes it `UNKNOWN`
  (flagged), not blank.
- **`test_serology_no_independent_cross_product`** — §6.6: the naive union
  cross-product `a/b/c + a/b/c` asserts excluded `a+a`/`b+b`; the engine must
  emit the exact set (here the valid factorization `a/c + b/c`, or exact `|`
  when a set does not factor), never the naive union.

---

## 9. What this shares with the sequence backend (summary)

| Concern | Sequence backend | Serology backend |
|---------|------------------|------------------|
| signature | base per interrogated column | projected specificity set |
| combined observation | per-column base **set** (diploid; §2.3.1) | union of specificity sets |
| candidate comparison | three-valued (main §4.1) | three-valued (same) |
| reference artifact | column model (D1) | serology model (§2) |
| missing-data symbol | `indeterminate` column (§4, D2) | `UNKNOWN` confidence (§5) |
| blank/absence | sequenced `gap` (D2) | `NOT_EXPRESSED` blank (§5, D7) |
| genotype ambiguity algo | shared (main §3.2–§3.3) | **same shared code** |
| `/` vs `|` reduction (exact default) | shared (main §3.3) | **same shared code** |
| phase regrouping | shared; D8 (product of molecules, not loci) | **same shared code**; D8 applies to the compatible-DNA set over phase-grouped loci |
| representation policy | D4 (`exact` default) | D4 |
| GL String content | DNA alleles | **DNA alleles only** (antigen report is separate, §7) |

The fact that a backend this different (specificities, broad/split lattice,
blanks) plugs into the same engine with only two methods overridden is the
evidence that the M2 abstraction (decision D6) is correct.

Nothing above is implemented; this is the agreed basis for the serology portion
of M2 once review signs off.
