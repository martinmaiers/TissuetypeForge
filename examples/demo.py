"""Working demonstration of TissueTypeForge Milestone 1.

Run with:  uv run python examples/demo.py

Shows, from one known truth GL String:
  * idealized first-field vs second-field observation (information loss)
  * loss of interlocus B~C phase under an unphased strategy
  * retention of B~C phase under a phase-aware strategy (spec §6)
  * a mixed-resolution registry strategy with a missing (untyped) locus
and prints the structured ObservationResult for one case.
"""

from __future__ import annotations

from pathlib import Path

from tissuetypeforge import load_strategy_file, observe, parse_gl_string, serialize

STRAT_DIR = Path(__file__).resolve().parents[1] / "strategies"

# Known phased diplotype (truth). B~C are phased in cis within each haplotype.
TRUTH = (
    "HLA-A*01:01:01:01+HLA-A*02:01:01^"
    "HLA-B*07:02:01~HLA-C*07:02:01+HLA-B*08:01:01~HLA-C*07:01:01^"
    "HLA-DQB1*02:01:01+HLA-DQB1*03:01:01"
)


def _rule() -> None:
    print("-" * 72)


def main() -> None:
    print("TissueTypeForge — Milestone 1 demonstration")
    _rule()
    print("TRUTH (known diplotype):")
    print(" ", TRUTH)
    # Confirm the parser/serializer round-trips the truth exactly.
    assert serialize(parse_gl_string(TRUTH)) == TRUTH
    print("  [GL String round-trips exactly]")
    _rule()

    for fname in (
        "first_field.toml",
        "second_field.toml",
        "second_field_bc_phased.toml",
        "mixed_registry.toml",
    ):
        strat = load_strategy_file(STRAT_DIR / fname)
        result = observe(TRUTH, strat, seed=1)
        print(f"strategy: {strat.name}")
        print(f"  observed : {result.observed}")
        print(f"  typed    : {', '.join(result.metadata.typed_loci) or '(none)'}")
        if result.metadata.missing_loci:
            print(f"  missing  : {', '.join(result.metadata.missing_loci)}")
        if result.metadata.retained_phase_groups:
            groups = "; ".join(
                "~".join(g) for g in result.metadata.retained_phase_groups
            )
            print(f"  phase    : retained {groups}")
        else:
            print("  phase    : no interlocus phase retained")
        print(f"  error    : {result.error_category.value}")
        _rule()

    # Full structured output for the phase-aware strategy.
    strat = load_strategy_file(STRAT_DIR / "second_field_bc_phased.toml")
    result = observe(TRUTH, strat, seed=1)
    print("Structured result (second_field_bc_phased):")
    print(f"  true_diplotype       : {result.true_diplotype}")
    print(f"  expected_observation : {result.expected_observation}")
    print(f"  observed             : {result.observed}")
    print(f"  error_category       : {result.error_category.value}")
    print(f"  reference_databases  : {result.metadata.reference_databases}")
    print("  information loss record:")
    for note in result.loss_record:
        print(f"    - {note}")


if __name__ == "__main__":
    main()
