"""Milestone 2 demonstration — strategy-file-driven observation.

Run with:  uv run python examples/demo_m2.py

Shows the full public-API path: a true GL String + a TOML strategy file +
a ReferenceRegistry of versioned artifacts -> observe() -> observed GL String,
a structured serologic antigen report, and information-loss metadata. No
internal engine functions are called.
"""

from __future__ import annotations

from pathlib import Path

from tissuetypeforge import load_strategy_file, observe
from tissuetypeforge.reference import Base, ColumnModel, FeatureKey, ReferenceRegistry
from tissuetypeforge.serology import (
    Assignment,
    BroadSplitLattice,
    Confidence,
    SerologyModel,
    Specificity,
)

STRAT = Path(__file__).resolve().parents[1] / "strategies"
EX2, EX3 = FeatureKey("exon", 2), FeatureKey("exon", 3)


def build_registry() -> ReferenceRegistry:
    """A small, pinned, hand-authored set of reference artifacts (toy data).

    In production these are built by a data-prep tool from a pinned
    IPD-IMGT/HLA release; here they are inline so the demo is self-contained.
    """
    reg = ReferenceRegistry()
    # genomic column models (exon2 cols 0-1, exon3 cols 2-3)
    reg.register_column_model("TOY-GEN", ColumnModel(
        locus="HLA-B", database_version="TOY-GEN/1.0",
        feature_of_column=(EX2, EX2, EX3, EX3),
        rows={
            "HLA-B*07:02": (Base.G, Base.A, Base.C, Base.C),
            "HLA-B*07:09": (Base.G, Base.A, Base.C, Base.C),  # sibling, same signature
            "HLA-B*08:01": (Base.T, Base.C, Base.G, Base.G),
        },
    ))
    reg.register_column_model("TOY-GEN", ColumnModel(
        locus="HLA-C", database_version="TOY-GEN/1.0",
        feature_of_column=(EX2, EX2, EX3, EX3),
        rows={
            "HLA-C*07:01": (Base.A, Base.A, Base.T, Base.T),
            "HLA-C*07:02": (Base.C, Base.C, Base.A, Base.A),
            "HLA-C*07:18": (Base.C, Base.C, Base.A, Base.A),  # sibling
        },
    ))
    # serology model
    B7, B15, B62 = Specificity("B7"), Specificity("B15"), Specificity("B62")
    reg.register_serology_model("TOY-SER", SerologyModel(
        locus="HLA-B", database_version="TOY-SER/1.0",
        assignment={
            "HLA-B*07:02": Assignment(frozenset({B7}), Confidence.ASSIGNED),
            "HLA-B*15:01": Assignment(frozenset({B62}), Confidence.ASSIGNED),
            "HLA-B*15:02": Assignment(frozenset({Specificity("B75")}), Confidence.ASSIGNED),
        },
        lattice=BroadSplitLattice(broad_of={B62: B15, Specificity("B75"): B15, B7: B7}),
    ))
    Cw7 = Specificity("Cw7")
    reg.register_serology_model("TOY-SER", SerologyModel(
        locus="HLA-C", database_version="TOY-SER/1.0",
        assignment={
            "HLA-C*07:01": Assignment(frozenset({Cw7}), Confidence.ASSIGNED),
            "HLA-C*07:02": Assignment(frozenset({Cw7}), Confidence.ASSIGNED),
        },
        lattice=BroadSplitLattice(broad_of={Cw7: Cw7}),
    ))
    return reg


def _show(title: str, truth: str, strategy_file: str, reg) -> None:
    print("-" * 72)
    print(f"strategy file : {strategy_file}")
    s = load_strategy_file(STRAT / strategy_file)
    r = observe(truth, s, registry=reg)
    print(f"truth         : {truth}")
    print(f"observed      : {r.observed}")
    if r.metadata.retained_phase_groups:
        groups = "; ".join("~".join(g) for g in r.metadata.retained_phase_groups)
        print(f"retained phase: {groups}")
    for report in r.antigen_reports:
        specs = sorted(sp.name for sp in report.specificities)
        print(f"antigen report: {report.locus} {report.resolution_level.value} "
              f"-> {specs}  (NOT in the GL String)")
    if r.metadata.reference_databases:
        print(f"reference dbs : {', '.join(r.metadata.reference_databases)}")


def main() -> None:
    print("TissueTypeForge — Milestone 2 demonstration (strategy-file driven)")
    reg = build_registry()

    _show(
        "SBT exon2+3 with B~C joint phase",
        "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01",
        "sbt_exon23_bc_phased.toml", reg,
    )
    _show(
        "Broad-only serology at B",
        "HLA-B*15:01+HLA-B*07:02",
        "serology_broad_B.toml", reg,
    )
    _show(
        "Mixed assay classes across loci",
        "HLA-A*01:01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01^"
        "HLA-C*07:01+HLA-C*07:02^HLA-DQB1*02:01+HLA-DQB1*03:01",
        "mixed_assay_registry.toml", reg,
    )
    print("-" * 72)


if __name__ == "__main__":
    main()
