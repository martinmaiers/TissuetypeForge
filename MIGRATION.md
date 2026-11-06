# Migration notes — 0.3.0

Version 0.3.0 deliberately **simplifies** the observation model and introduces
population-level joint typing-profile distributions. This is an intentional
breaking change: the model was previously more complex than our initial
validation experiments require.

## Summary of intentional breaking changes

| Area | Before (0.1.x / 0.2.x) | After (0.3.0) |
|------|------------------------|---------------|
| Resolution vocabulary | `allele`, `first_field`, `second_field`, `serology`, `intermediate`, `full_gene` | exactly `0-field`, `1-field`, `2-field` |
| Config unit | `Strategy` (per-locus `assay_class`, resolution, reference, phase) | `TypingProfile` + `PopulationTypingModel` |
| Per-locus sampling | independent per-locus resolution | **joint** — one complete profile sampled per donor |
| Interlocus phase | `phase_groups` (e.g. B~C) retained | removed; observations never assert interlocus phase |
| Assays | serology / SBT / SSO / idealized, reference registries | none — idealized field truncation only |
| Entry point | `observe(truth_gl, strategy, *, registry=...)` | `observe_genotype(truth_gl, profile)` and `simulate_cohort(...)` |
| Result type | `ObservationResult` + `ObservationMetadata` | `DonorObservation`, `CohortResult`, `ReproducibilityManifest` |

### Removed modules (superseded)

- `tissuetypeforge.observation` (framework, sequence_engine, joint_phase, policy,
  reduce, resolution)
- `tissuetypeforge.serology`
- `tissuetypeforge.reference`
- `tissuetypeforge.strategy`
- `tissuetypeforge.result`

### Retained general-purpose functionality

- `tissuetypeforge.glstring` — the GL String wrapper over the SIP `glstring`
  library (`parse_gl_string`, `serialize`, AST node types). Unchanged.
- `tissuetypeforge.nomenclature.truncate_to_fields` — field truncation helper.

## Porting a previous configuration

### Old: per-locus `Strategy` with B~C phasing (TOML)

```toml
name = "second_field_bc_phased"
reference_database = "IPD-IMGT/HLA 3.55.0"

[defaults]
resolution = "second_field"

[loci.HLA-C]
resolution = "first_field"

[[phase_groups]]
loci = ["B", "C"]
```

### New: a typing profile inside a population (TOML)

Map the old resolution modes to the new vocabulary:

| Old mode | New mode |
|----------|----------|
| `first_field` | `1-field` |
| `second_field` | `2-field` |
| a locus that was `typed = false` | `0-field` |
| `allele` / `serology` / `intermediate` / `full_gene` | *not supported* — choose `2-field` (closest retained resolution) or model as out of scope |

The `phase_groups` declaration is **dropped** — interlocus phase is no longer
modeled.

```toml
name = "registry-example"
loci = ["A", "C", "B", "DRB1"]

# The old single strategy becomes ONE profile. Add more profiles (each with a
# weight) to model a population typed under multiple historical protocols.
[[profiles]]
id = "second_field_bc"   # phase is not retained; the id is just a label
weight = 1.0
resolution = { A = "2-field", B = "2-field", C = "1-field", DRB1 = "2-field" }
```

### Old call site

```python
from tissuetypeforge import load_strategy_file, observe

strat = load_strategy_file("strategies/second_field.toml")
result = observe(truth_gl, strat, seed=1)
result.observed
```

### New call site

```python
from tissuetypeforge import load_population_file, observe_genotype, simulate_cohort

pop = load_population_file("populations/registry_example.toml")

# One donor under one profile:
obs = observe_genotype(truth_gl, pop.profiles[0])
obs.observed

# A whole cohort with reproducible, seeded profile sampling:
result = simulate_cohort([truth_gl, ...], pop, seed=1)
result.manifest.profile_assignments
```

## Why these changes

- **Three resolution modes** cover the information-loss behaviors our validation
  experiments need, without the maintenance burden of assay-physics modeling.
- **Joint profiles** reflect reality: a donor's loci are typed by one protocol,
  so their resolutions are correlated. Independent per-locus sampling produced
  combinations no real lab would report.
- **No interlocus phase**: the initial experiments do not consume phase, and
  retaining it added substantial complexity (`joint_phase`, phase groups).
- **No serology/SBT/SSO**: explicitly out of scope for this milestone.

If you need the retired assay-derived ambiguity, serology, or phasing behavior,
pin to a `0.2.x` release.
