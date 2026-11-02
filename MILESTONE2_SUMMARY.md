# Milestone 2 Summary — Assay-Derived Ambiguity and Serology

**Status: ready for acceptance.** 88 tests pass, `ruff check .` clean, both
demos run. No stochastic mistyping (M3) or batch-scale simulation (M4) was
built. No remote repositories, commits, or published packages were created.

## Scope delivered

M2 adds scientifically realistic, assay-derived ambiguity on top of the M1
observation spine, driven entirely by the public API (a truth GL String + a
TOML strategy + a reference registry → `observe()`), with no internal engine
calls required of the user.

| Area | Delivered |
|------|-----------|
| GL String | Migrated to the published SIP `glstring` library (GL String 1.1); internal parser retired |
| Reference artifacts (M2A) | Versioned, checksummed per-locus `ColumnModel`; `gap` vs `indeterminate` as distinct symbols; serology `SerologyModel` with broad/split lattice |
| Sequence engine (M2B) | Diploid **pair enumeration**, three-valued candidate comparison, `reference_policy`, exact `\|` / lossless `/`, resource limits, whole-molecule joint phase |
| Serology backend (M2C) | Separate antigen report + compatible-DNA GL String; broad/split; declared Bw4/Bw6 epitopes; blank vs unknown |
| Strategy integration (D4) | `assay_class` dispatch (idealized / sbt / sso / serology), reference selection by name, phase-observation groups, strict/permissive policy, resource limits, mixed assays per locus |

## Final decisions, as implemented

- **D1 Public epitopes** — `observe_serology(..., observed_epitopes=...)` and the
  strategy's per-locus `observed_epitopes` list. Only explicitly declared
  epitopes participate; unobserved Bw4/Bw6 never narrows the candidate set. The
  reference-model assignment (all epitopes an allele carries) is kept distinct
  from the simulated observation (the declared subset).
- **D2 GL factorization** — bounded `3^n` bipartition search for lossless `/`;
  exact `|` enumeration otherwise. Resource limits raise explicit
  `ResourceLimitError`. No further optimization spent.
- **D3 Joint-phase scope** — two-molecule diploid model; member loci may span
  multiple interrogated regions (B~C over exon 2 **and** exon 3 is tested).
  Verified: candidate haplotypes reproduce the joint observation; phase retained
  only where observed; no fabricated cross-locus haplotypes; multi-region
  ambiguity handled consistently.
- **D4 Strategy integration** — `reference_policy` defaults to `strict`
  (fail-closed via `IndeterminateReferenceError`); `permissive` flags
  approximate results; `sequence_supplied` does not erase candidate-side
  uncertainty. Resource limits are configurable per strategy.

## Scientific guarantees (enforced by tests)

1. **Diploid genotype ambiguity, not independent matching.** Candidate pairs are
   enumerated against the combined diploid observation; a `/`-cross-product that
   would assert excluded genotypes is never emitted (sequence §7.2, serology
   §6.6).
2. **Reference knowledge ≠ laboratory signal.** Incomplete reference → an
   *indeterminate* comparison. `strict` raises rather than silently excluding;
   `permissive` keeps and flags. A sequenced `gap` is discriminating; missing
   data is not.
3. **Phase only where observed.** A phase-resolving `[[phase_groups]]` enumerates
   whole-molecule pairings (product of molecules, not of loci). An unphased
   strategy introduces no cross-locus phase.
4. **Serology is two outputs.** Antigen/epitope tokens live in a structured
   report and are never placed in a GL String; the compatible-DNA GL String
   carries only nomenclature alleles.
5. **Backward compatible.** An all-idealized strategy needs no registry and
   produces identical output to M1.

## Dependency on the published SIP library

All GL String parsing/serialization uses the published **`glstring`** package
(PyPI `glstring`, repo `societyforimmunepolymorphism/py-glstring`), a
standards-conformant GL String 1.1 implementation. TissueTypeForge no longer
ships a competing parser; `src/tissuetypeforge/glstring/__init__.py` is a thin
adapter. Declared in `pyproject.toml` as `glstring>=0.1.0`.

## Public API

```python
from tissuetypeforge import (
    observe, load_strategy_file, load_strategy,
    Strategy, ObservationResult, ErrorCategory,
    parse_gl_string, serialize,
)
from tissuetypeforge.reference import ReferenceRegistry, ColumnModel, Base, FeatureKey
from tissuetypeforge.serology import (
    SerologyModel, Assignment, Confidence, Specificity, Epitope,
)
```

`observe(truth_gl, strategy, *, seed=None, registry=None) -> ObservationResult`
with fields `observed` (GL String, DNA only), `antigen_reports` (structured),
`metadata` (typed/missing loci, retained phase, reference versions, per-locus
notes), `expected_observation`, `error_category`, `loss_record`.

## End-to-end acceptance tests (the 12 required demonstrations)

All in `tests/test_acceptance.py`, driven through the public API:

1. Idealized M1 typing unchanged — `test_1_idealized_unchanged`
2. Sequence ambiguity via TOML — `test_2_sequence_ambiguity_via_toml`
3. Diploid cis/trans exact alternatives — `test_3_cis_trans_exact_alternatives`
4. Serology report + separate DNA GL — `test_4_serology_report_and_dna_gl`
5. Broad-only ignores unobserved Bw4/Bw6 — `test_5_broad_only_ignores_unobserved_epitopes`
6. Declared epitopes constrain genotypes — `test_6_declared_epitopes_constrain`
7. B~C joint phase → 4 diplotypes — `test_7_bc_joint_phase_four_diplotypes`
   (+ `test_7b_joint_phase_multi_region_no_fabricated_haplotype`, D3)
8. Unphased → no cross-locus phase — `test_8_unphased_no_cross_locus_phase`
9. Mixed-assay strategy across loci — `test_9_mixed_assay_strategy`
10. Strict reference fails on incompleteness — `test_10_strict_fails_on_incomplete_reference`
11. Error-free observation compatible with truth — `test_11_error_free_compatible_with_truth`
12. `glstring` round-trips all outputs — `test_12_glstring_roundtrips_all_outputs`

## Verification

```
uv run pytest        # 88 passed
uv run ruff check .  # All checks passed!
uv run python examples/demo.py       # M1
uv run python examples/demo_m2.py    # M2 (strategy-file driven)
```

Worked-example arithmetic (§7.2 sequence cis/trans, §6.4.1 joint phase 4-vs-8,
§6.2 broad-split four-way, §6.6 serology cross-product) was independently
verified numerically during development and is asserted in the suite.

## Acceptance recommendation

M2 meets the stated scope and all four decisions, every required acceptance test
passes, the suite and linter are clean, and the public API lets a user go from a
true GL String + strategy file to an observed typing without touching internal
engines. **M2 is ready for acceptance.**

### Noted for M3 / future review

- The joint-phase engine is the diploid two-molecule model (as decided);
  multi-region per locus is supported and tested, but >2-copy systems (KIR
  copy-number) remain future work.
- `resolves_dosage` is reserved in the capability vocabulary but inert in M2
  (clean diploid); it becomes meaningful only with the M3 mixture model.
- Error models (dropout, false homozygosity, substitution) and batch-scale
  simulation are deferred to M3/M4 as instructed.
