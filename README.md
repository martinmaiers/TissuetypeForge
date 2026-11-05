# TissueTypeForge

Simulate HLA typing observations from known genetic truth.

TissueTypeForge is an open-source Python package of the **Society for Immune
Polymorphism (SIP)**, part of the **HLA Forge** initiative. Its companion
project, **DiplotypeForge**, generates synthetic phased diplotypes from
haplotype frequencies. TissueTypeForge is independently usable: it accepts real
or simulated known genotypes without any dependency on DiplotypeForge.

> **Core principle:** simulate what the typing method *observes*, including what
> it *cannot* observe.

Given a known HLA genotype (a GL String) and a **typing profile**,
TissueTypeForge simulates the typing result that profile would report —
reproducing resolution-driven information loss and missing results. It also
models a **population** as a joint distribution of typing profiles, reflecting
how real donor registries contain groups of donors typed under different
historical protocols.

## The model

### 1. Controlled resolution vocabulary

Exactly three resolution modes are supported, applied per locus:

| Mode | Meaning | `A*02:01:01:01 + A*03:01:01:01` becomes |
|------|---------|------------------------------------------|
| `0-field` | no typing available at this locus | *(locus omitted — a missing result)* |
| `1-field` | retain only the first field | `HLA-A*02+HLA-A*03` |
| `2-field` | retain the first two fields | `HLA-A*02:01+HLA-A*03:01` |

Heterozygosity is preserved where distinguishable at the requested resolution.
If two different alleles collapse to the same first-field name, the result does
**not** assert true homozygosity: the collapse is recorded as preserved
uncertainty in the observation's `loss_record`.

Serology, SBT, SSO, SSP, and intermediate-resolution methods are **not** part of
this model.

### 2. Typing profiles and populations

A **typing profile** assigns a resolution mode to *every* configured locus and
carries a unique id and a weight. A **population typing model** is a categorical
distribution over profiles. For each simulated donor, exactly one profile is
sampled and applied as a whole — so resolution choices across loci are **joint**,
never sampled independently per locus. This models correlated, protocol-driven
typing patterns (e.g. a "legacy" lab types every locus coarsely).

The population's *genetic composition* (which true genotypes occur) and its
*typing-profile distribution* (how donors were typed) are kept as separate
concepts. No association between genotype and typing profile is introduced.

### 3. No interlocus phase

Observations never assert interlocus phase. Each typed locus is reported as an
independent `^`-delimited block. Intralocus allele ambiguity (`/`) is retained.

## Public API

```python
from tissuetypeforge import (
    # resolution vocabulary
    ResolutionMode, parse_resolution,
    # profiles / population
    TypingProfile, PopulationTypingModel,
    build_profile, build_population,
    load_population, load_population_file,
    # single-donor observation
    observe_genotype, DonorObservation,
    # cohort simulation
    simulate_cohort, sample_profile,
    CohortResult, ReproducibilityManifest,
    # general-purpose GL String round-trip
    parse_gl_string, serialize,
)
```

### Observe one donor under one profile

```python
from tissuetypeforge import TypingProfile, ResolutionMode, observe_genotype

profile = TypingProfile(
    id="modern",
    weight=1.0,
    resolution={"A": ResolutionMode.TWO_FIELD, "C": ResolutionMode.TWO_FIELD},
)
obs = observe_genotype("HLA-A*02:01:01:01+HLA-A*03:01:01:01", profile)
obs.observed        # "HLA-A*02:01+HLA-A*03:01"
obs.profile_id      # "modern"
obs.missing_loci    # ()
```

### Simulate a cohort from a population config

```toml
# populations/registry_example.toml
name = "registry-example"
loci = ["A", "C", "B", "DRB1"]

[[profiles]]
id = "legacy"
weight = 0.20
resolution = { A = "1-field", B = "1-field", C = "0-field", DRB1 = "1-field" }

[[profiles]]
id = "intermediate"
weight = 0.50
resolution = { A = "2-field", B = "2-field", C = "1-field", DRB1 = "2-field" }

[[profiles]]
id = "modern"
weight = 0.30
resolution = { A = "2-field", B = "2-field", C = "2-field", DRB1 = "2-field" }
```

```python
from tissuetypeforge import load_population_file, simulate_cohort

pop = load_population_file("populations/registry_example.toml")
result = simulate_cohort(truth_genotypes, pop, seed=42, donor_ids=donor_ids)

for obs in result.observations:
    print(obs.profile_id, obs.observed)

result.manifest.assigned_profile_counts   # profile id -> donor count
result.manifest.profile_assignments       # per-donor profile ids (replayable)
result.manifest.seed                       # the seed used
```

The same inputs, configuration, and seed always produce identical outputs.

### Batch driver (file of `id,glstring`)

A command-line driver processes a file of donor records — a donor id and that
donor's *true* genotype as a GL String — and writes the simulated observed
typings plus a reproducibility manifest:

```bash
# examples/donors.csv has an `id,glstring` header
uv run python -m tissuetypeforge.driver examples/donors.csv \
    --population populations/registry_example.toml \
    --seed 42 --out observed.tsv --manifest observed.manifest.json
```

The output TSV has columns `donor_id, profile_id, true_genotype, observed,
missing_loci, status, error_detail`. Malformed records are isolated at the
record level (`--mode permissive`, the default) or abort the batch
(`--mode strict`). A JSON manifest records the seed, population, profile
weights, per-donor profile assignments, and input/output SHA-256 checksums. The
installed package also exposes this as the `tissuetypeforge` console script.

## Configuration validation

`build_population` / `load_population(_file)` validate before simulation and
raise `ProfileValidationError` with a clear message when:

- profile ids are not unique;
- a weight is not finite or is negative;
- no profile has positive weight;
- weights do not sum to 1 (within a small tolerance);
- a profile does not specify every configured locus;
- a resolution value is outside the `0-field` / `1-field` / `2-field` vocabulary;
- a profile names an unknown locus.

## GL String dependency

TissueTypeForge uses the published **SIP `glstring`** library (PyPI `glstring`)
for all GL String parsing and serialization — a standards-conformant GL String
1.1 implementation. TissueTypeForge does **not** ship its own GL String parser.

## Development

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run python examples/population_demo.py
```

## Migration from earlier versions

The observation model was deliberately simplified in 0.3.0 (new resolution
vocabulary, population-level joint profiles, no interlocus phasing, no
serology/SBT/SSO). See [MIGRATION.md](MIGRATION.md) for the breaking changes and
how to port previous `Strategy`/`phase_group` configurations.

## License

Apache-2.0 (see `LICENSE`).
