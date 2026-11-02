# TissueTypeForge

Simulate HLA typing observations from known genetic truth.

TissueTypeForge is an open-source Python package of the **Society for Immune
Polymorphism (SIP)**, part of the **HLA Forge** initiative. Its companion
project, **DiplotypeForge**, generates synthetic phased diplotypes from
haplotype frequencies. TissueTypeForge is independently usable: it accepts real
or simulated known genotypes without any dependency on DiplotypeForge.

> **Core principle:** simulate what the typing method *observes*, including what
> it *cannot* observe.

Given a known phased HLA diplotype (a GL String) and a configurable typing
strategy, TissueTypeForge simulates the typing result that strategy would
report — reproducing information loss, ambiguity, and (in later milestones)
mistyping.

## Status: Milestone 1 — Observation framework

This release implements the observation *spine* on which all later assay models
are built. It is deliberately scoped:

| Area | M1 scope |
|------|----------|
| GL String | Standards-conformant recursive-descent parser + serializer, round-trippable |
| Internal model | Phase-aware genotype representation supporting arbitrary locus subsets with known phase |
| Strategy schema | Machine-readable TOML, per-locus config, phase groups, validation |
| Resolution | **Idealized** `allele`, `first_field`, `second_field` observations |
| Coverage | Per-locus typed / not-typed (missing results) |
| Phase | Explicit phase-retention rules mapping to GL String operators |
| Output | Observed GL String + structured metadata + information-loss record |

### Explicitly deferred (not in M1)

To avoid presenting idealized transformations as empirical laboratory
simulations, the following raise an explicit `NotImplementedError` rather than
being faked:

- Sequence-region / assay-derived ambiguity (**M2**)
- Serology mappings (**M2**)
- Intermediate-resolution models (**M2**)
- Stochastic error models: dropout, false homozygosity, substitution (**M3**)
- Registry-scale batch operation (**M4**)

The `ObservationResult` already exposes `expected_observation`, `observed`, and
`error_category` fields so the M3 error layer slots in without changing the
output contract. In M1, `expected_observation == observed` and the error
category is always `none`.

## Scientific assumptions (M1)

- **Idealized first/second-field** observation is modeled as nomenclature-field
  truncation of the *true* allele. This is a baseline convenience model, **not**
  an empirical assay. It assumes the reported field-level name is exactly the
  truncation of truth with no measurement ambiguity. Real assays interrogate
  sequence and are modeled in M2.
- **Phase is only retained when a strategy declares a phase group** containing
  the relevant loci. Phase is never invented because truth is known
  (spec §6). Loci not sharing a phase group are reported with no interlocus
  phase asserted.
- **Idealized modes do not introduce ambiguity.** A field-truncated allele is
  reported as a single allele, not a MAC/ambiguity set, because the idealized
  model asserts the method resolves exactly to that field. Assay-derived
  ambiguity is an M2 concern.

## GL String operator mapping

TissueTypeForge uses the GL String operators with standard precedence
(tightest to loosest):

| Operator | Meaning | Internal concept |
|----------|---------|------------------|
| `/` | allele ambiguity | `AlleleList` |
| `~` | in-phase (cis / haplotype) | ordering within a `GenotypeBlock` |
| `+` | gene-copy pairing (trans) | the two copies of a `GenotypeBlock` |
| `\|` | genotype ambiguity | `GenotypeAmbiguity` |
| `^` | locus delimiter | separate `GenotypeBlock`s |

A **phased block** over loci `{B, C}` is a single `GenotypeBlock` serialized as
`B*..~C*..+B*..~C*..`. Breaking phase between B and C yields two independent
blocks `B*..+B*.. ^ C*..+C*..`, which asserts no interlocus relationship.

## Milestone 2 — assay-derived ambiguity and serology

M2 adds scientifically realistic, assay-derived ambiguity on top of the M1
spine. A strategy now selects an **assay class** per locus and the engine
simulates what that assay can and cannot distinguish.

| Assay class | Behavior |
|-------------|----------|
| `idealized` | M1 nomenclature truncation (unchanged, no reference data) |
| `sbt` / `sso` | sequence observation over interrogated features; genotype-level ambiguity computed by **diploid pair enumeration** |
| `serology` | antigen observation via a versioned allele→specificity map, broad/split lattice, and explicitly declared public epitopes |

Key scientific guarantees:

- **Diploid genotype ambiguity, not independent allele matching.** The engine
  enumerates unordered allele *pairs* that reproduce the combined diploid
  observation, so it never fabricates a genotype by crossing per-allele matches.
- **Reference knowledge ≠ laboratory signal.** Incomplete reference data becomes
  an *indeterminate* comparison, never a silent match or a silent exclusion.
  `reference_policy = strict` (default) raises rather than returning a reduced
  set it cannot prove exact; `permissive` keeps indeterminate candidates and
  flags the result as approximate.
- **Phase only where observed.** A `[[phase_groups]]` of sequence-typed loci is
  enumerated as whole-molecule haplotype pairings (the product of observed
  molecules, not of per-locus ambiguities), so partial phase is reproduced
  without inventing cross-locus haplotypes.
- **Serology is two outputs.** A structured antigen report (broad/split
  specificities + declared Bw4/Bw6 epitopes) is kept **separate** from the
  compatible-DNA GL String; antigen tokens are never placed in a GL String.

### Public API

```python
from tissuetypeforge import (
    observe,                 # observe(truth_gl, strategy, *, seed=None, registry=None)
    load_strategy_file,      # load a Strategy from a TOML file
    load_strategy,           # load a Strategy from a parsed mapping
    Strategy, ObservationResult, ErrorCategory,
    parse_gl_string, serialize,   # thin wrappers over the glstring library
)
from tissuetypeforge.reference import ReferenceRegistry, ColumnModel, Base, FeatureKey
from tissuetypeforge.serology import SerologyModel, Assignment, Confidence, Specificity, Epitope
```

End-to-end usage (truth GL String + strategy file + reference registry):

```python
registry = ReferenceRegistry()
registry.register_column_model("IPD-3.55", column_model_for_B)   # pinned artifact
strategy = load_strategy_file("strategies/sbt_exon23_bc_phased.toml")
result = observe("HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01",
                 strategy, registry=registry)
result.observed          # observed GL String (DNA alleles only)
result.antigen_reports   # structured serologic antigen reports (if any)
result.metadata          # typed/missing loci, retained phase, reference versions
```

See `examples/demo_m2.py` for a runnable, strategy-file-driven demonstration of
SBT joint phase, broad-only serology, and a mixed-assay strategy.

### GL String dependency

TissueTypeForge uses the published **SIP `glstring`** library
(repo `societyforimmunepolymorphism/py-glstring`, PyPI `glstring`) for all GL
String parsing and serialization — a standards-conformant GL String 1.1
implementation. TissueTypeForge does **not** ship its own GL String parser. The
dependency is declared in `pyproject.toml` (`glstring>=0.1.0`).

## Development

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run python examples/demo.py       # M1
uv run python examples/demo_m2.py    # M2
```

## License

Apache-2.0 (see `LICENSE`).
