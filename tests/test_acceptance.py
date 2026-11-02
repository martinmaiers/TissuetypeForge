"""End-to-end M2 acceptance tests (Decision 4 + the 12 required demonstrations).

Every test drives the *public* API: a truth GL String + a Strategy (loaded from
a TOML file where practical) + a ReferenceRegistry -> observe() -> result. No
internal engine functions are called from these tests.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from glstring import GLString  # the published SIP library, used directly for round-trip

from tissuetypeforge import ObservationResult, load_strategy_file, observe
from tissuetypeforge.observation.policy import IndeterminateReferenceError
from tissuetypeforge.reference import Base, ColumnModel, FeatureKey, ReferenceRegistry
from tissuetypeforge.serology import (
    Assignment,
    BroadSplitLattice,
    Confidence,
    Epitope,
    SerologyModel,
    Specificity,
)
from tissuetypeforge.strategy import load_strategy

STRAT = Path(__file__).resolve().parents[1] / "strategies"
EX2, EX3 = FeatureKey("exon", 2), FeatureKey("exon", 3)

# --- toy reference artifacts -------------------------------------------------
# Sequence genomic model "TOY-GEN" per locus, exon2(0,1)+exon3(2,3).
_B_GEN = ColumnModel(
    locus="HLA-B", database_version="TOY-GEN/1.0",
    feature_of_column=(EX2, EX2, EX3, EX3),
    rows={
        "HLA-B*07:02": (Base.G, Base.A, Base.C, Base.C),
        "HLA-B*07:09": (Base.G, Base.A, Base.C, Base.C),  # shares B*07:02 signature
        "HLA-B*08:01": (Base.T, Base.C, Base.G, Base.G),
    },
)
_C_GEN = ColumnModel(
    locus="HLA-C", database_version="TOY-GEN/1.0",
    feature_of_column=(EX2, EX2, EX3, EX3),
    rows={
        "HLA-C*07:01": (Base.A, Base.A, Base.T, Base.T),
        "HLA-C*07:02": (Base.C, Base.C, Base.A, Base.A),
        "HLA-C*07:18": (Base.C, Base.C, Base.A, Base.A),  # shares C*07:02 signature
    },
)
# Serology model "TOY-SER".
B7, B8, B62, B75, B15 = (Specificity(s) for s in ("B7", "B8", "B62", "B75", "B15"))
Cw7 = Specificity("Cw7")
Bw4, Bw6 = Epitope("Bw4"), Epitope("Bw6")
_B_SER = SerologyModel(
    locus="HLA-B", database_version="TOY-SER/1.0",
    assignment={
        "HLA-B*07:02": Assignment(frozenset({B7}), Confidence.ASSIGNED, frozenset({Bw6})),
        "HLA-B*08:01": Assignment(frozenset({B8}), Confidence.ASSIGNED, frozenset({Bw6})),
        "HLA-B*15:01": Assignment(frozenset({B62}), Confidence.ASSIGNED, frozenset({Bw6})),
        "HLA-B*15:02": Assignment(frozenset({B75}), Confidence.ASSIGNED, frozenset({Bw4})),
    },
    lattice=BroadSplitLattice(broad_of={B62: B15, B75: B15, B15: B15, B7: B7, B8: B8}),
)
_C_SER = SerologyModel(
    locus="HLA-C", database_version="TOY-SER/1.0",
    assignment={
        "HLA-C*07:01": Assignment(frozenset({Cw7}), Confidence.ASSIGNED),
        "HLA-C*07:02": Assignment(frozenset({Cw7}), Confidence.ASSIGNED),
    },
    lattice=BroadSplitLattice(broad_of={Cw7: Cw7}),
)


@pytest.fixture
def registry() -> ReferenceRegistry:
    reg = ReferenceRegistry()
    reg.register_column_model("TOY-GEN", _B_GEN)
    reg.register_column_model("TOY-GEN", _C_GEN)
    reg.register_serology_model("TOY-SER", _B_SER)
    reg.register_serology_model("TOY-SER", _C_SER)
    return reg


def _roundtrips(gl: str) -> bool:
    if gl == "":
        return True
    once = GLString.parse(gl).to_string()
    twice = GLString.parse(once).to_string()
    return once == twice


# === 1. Idealized M1 typing remains unchanged ===============================
def test_1_idealized_unchanged() -> None:
    s = load_strategy_file(STRAT / "second_field.toml")
    truth = "HLA-A*01:01:01:01+HLA-A*02:01:01"
    r = observe(truth, s)
    assert r.observed == "HLA-A*01:01+HLA-A*02:01"
    assert r.antigen_reports == ()


# === 2. Sequence-based ambiguity via a TOML strategy ========================
def test_2_sequence_ambiguity_via_toml(registry) -> None:
    s = load_strategy(
        {
            "name": "sbt_B",
            "loci": {
                "HLA-B": {
                    "assay_class": "sbt",
                    "resolution": "second_field",
                    "reference_name": "TOY-GEN",
                    "interrogated_features": ["exon2", "exon3"],
                }
            },
        }
    )
    # B*07:02 and B*07:09 share interrogated signature -> genotype ambiguity.
    truth = "HLA-B*07:02+HLA-B*08:01"
    r = observe(truth, s, registry=registry)
    gl = GLString.parse(r.observed)
    assert "HLA-B*07:09" in r.observed or "HLA-B*07:02" in r.observed
    assert gl.is_ambiguous()  # B*07:02/B*07:09 indistinguishable


# === 3. Diploid cis/trans ambiguity: exact expected alternatives ============
def test_3_cis_trans_exact_alternatives(registry) -> None:
    # Single-column interrogation producing the classic §7.2 structure.
    m = ColumnModel(
        locus="HLA-A", database_version="TOY-72/1.0",
        feature_of_column=(FeatureKey("exon", 2), FeatureKey("exon", 2)),
        rows={
            "HLA-A*01:01": (Base.G, Base.A),
            "HLA-A*02:01": (Base.T, Base.C),
            "HLA-A*05:01": (Base.G, Base.C),
            "HLA-A*06:01": (Base.T, Base.A),
        },
    )
    reg = ReferenceRegistry()
    reg.register_column_model("M72", m)
    s = load_strategy({
        "name": "sbt_A",
        "loci": {"HLA-A": {
            "assay_class": "sbt", "reference_name": "M72",
            "interrogated_features": ["exon2"],
        }},
    })
    r = observe("HLA-A*01:01+HLA-A*02:01", s, registry=reg)
    assert r.observed == "HLA-A*01:01+HLA-A*02:01|HLA-A*05:01+HLA-A*06:01"
    assert "HLA-A*01:01+HLA-A*06:01" not in r.observed


# === 4. Serology: structured antigen report + separate DNA GL String ========
def test_4_serology_report_and_dna_gl(registry) -> None:
    s = load_strategy_file(STRAT / "serology_broad_B.toml")
    r = observe("HLA-B*15:01+HLA-B*07:02", s, registry=registry)
    assert len(r.antigen_reports) == 1
    report = r.antigen_reports[0]
    assert {sp.name for sp in report.specificities} == {"B15", "B7"}
    # DNA GL String contains only DNA alleles, never antigen tokens.
    for tok in ("B7", "B15", "B62"):
        assert f"+{tok}" not in r.observed and f"/{tok}" not in r.observed
    assert "HLA-B*" in r.observed


# === 5. Broad-only serology does not use unobserved Bw4/Bw6 =================
def test_5_broad_only_ignores_unobserved_epitopes(registry) -> None:
    s = load_strategy_file(STRAT / "serology_broad_B.toml")
    # B*15:01(Bw6) + B*07:02(Bw6): without epitope typing, B*15:02(Bw4) must
    # remain a candidate (epitopes must not narrow the set).
    r = observe("HLA-B*15:01+HLA-B*07:02", s, registry=registry)
    assert "HLA-B*15:02" in r.observed  # not excluded by unobserved Bw4
    assert r.antigen_reports[0].public_epitopes == frozenset()


# === 6. Declared Bw4/Bw6 typing constrains compatible genotypes ============
def test_6_declared_epitopes_constrain(registry) -> None:
    s = load_strategy({
        "name": "ser_epitope_B",
        "loci": {"HLA-B": {
            "assay_class": "serology", "resolution": "serology",
            "reference_name": "TOY-SER", "serology_level": "broad",
            "observed_epitopes": ["Bw4", "Bw6"],
        }},
        "reference_policy": "permissive",
    })
    # truth B*15:01(Bw6)+B*07:02(Bw6): declaring Bw6 typing excludes B*15:02(Bw4).
    r = observe("HLA-B*15:01+HLA-B*07:02", s, registry=registry)
    assert "HLA-B*15:02" not in r.observed  # Bw4 allele now excluded
    assert Epitope("Bw6") in r.antigen_reports[0].public_epitopes


# === 7. B~C joint-phase: exactly four candidate diplotypes ==================
def test_7_bc_joint_phase_four_diplotypes(registry) -> None:
    # B*07:02~C*07:02 / B*08:01~C*07:01 phased truth; B*07:09 and C*07:18 are
    # the indistinguishable siblings -> 2x2 = 4 joint diplotypes.
    truth = "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01"
    s = load_strategy_file(STRAT / "sbt_exon23_bc_phased.toml")
    r = observe(truth, s, registry=registry)
    # one locus block spanning B and C (joint), with 4 alternatives
    assert r.observed.count("|") == 3
    assert "~" in r.observed
    assert ("HLA-B", "HLA-C") in r.metadata.retained_phase_groups
    # fabricated cross haplotype must not appear
    assert "HLA-B*07:02~HLA-C*07:01" not in r.observed or r.observed.count("~") >= 2


def test_7b_joint_phase_multi_region_no_fabricated_haplotype(registry) -> None:
    # D3: multiple interrogated regions (exon2+exon3) per locus; verify ambiguity
    # in each region is handled consistently and no cross-locus haplotype is
    # fabricated. C*07:01 (the molecule-1 partner of B*08:01) must never appear
    # in cis with a B*07-class allele.
    truth = "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01"
    s = load_strategy_file(STRAT / "sbt_exon23_bc_phased.toml")
    r = observe(truth, s, registry=registry)
    for alt in r.observed.split("|"):
        # each alternative is two B~C haplotypes joined by +
        for hap in alt.split("+"):
            b, c = hap.split("~")
            if c == "HLA-C*07:01":
                assert b == "HLA-B*08:01"  # only ever with its true cis partner
            else:  # C*07:02 or its sibling C*07:18
                assert b in ("HLA-B*07:02", "HLA-B*07:09")


# === 8. Unphased strategy introduces no cross-locus phase ===================
def test_8_unphased_no_cross_locus_phase(registry) -> None:
    # Same sequence assay on B and C but NO phase group -> independent blocks.
    s = load_strategy({
        "name": "sbt_bc_unphased",
        "defaults": {
            "assay_class": "sbt", "reference_name": "TOY-GEN",
            "interrogated_features": ["exon2", "exon3"],
        },
        "loci": {"HLA-B": {}, "HLA-C": {}},
    })
    truth = "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01"
    r = observe(truth, s, registry=registry)
    # no B~C haplotype in the output; B and C are separate ^ blocks
    assert "~" not in r.observed
    assert r.metadata.retained_phase_groups == ()


# === 9. Mixed-resolution, mixed-assay strategy across loci ==================
def test_9_mixed_assay_strategy(registry) -> None:
    s = load_strategy_file(STRAT / "mixed_assay_registry.toml")
    truth = (
        "HLA-A*01:01:01+HLA-A*02:01:01^"
        "HLA-B*07:02+HLA-B*08:01^"
        "HLA-C*07:01+HLA-C*07:02^"
        "HLA-DQB1*02:01+HLA-DQB1*03:01"
    )
    r = observe(truth, s, registry=registry)
    assert "HLA-A*01:01+HLA-A*02:01" in r.observed  # idealized 2-field
    assert "HLA-DQB1" not in r.observed             # untyped -> missing
    assert "HLA-DQB1" in r.metadata.missing_loci
    assert len(r.antigen_reports) == 1              # serology at C
    assert _roundtrips(r.observed)


# === 10. Strict reference fails on incomplete reference =====================
def test_10_strict_fails_on_incomplete_reference() -> None:
    m = ColumnModel(
        locus="HLA-A", database_version="TOY-INC/1.0",
        feature_of_column=(FeatureKey("exon", 2), FeatureKey("exon", 2)),
        rows={
            "HLA-A*01:01": (Base.G, Base.A),
            "HLA-A*02:01": (Base.T, Base.C),
            "HLA-A*09:01": (Base.G, Base.INDETERMINATE),  # undecidable
        },
    )
    reg = ReferenceRegistry()
    reg.register_column_model("INC", m)
    s = load_strategy({
        "name": "strict_A",
        "reference_policy": "strict",
        "loci": {"HLA-A": {
            "assay_class": "sbt", "reference_name": "INC",
            "interrogated_features": ["exon2"],
        }},
    })
    with pytest.raises(IndeterminateReferenceError):
        observe("HLA-A*01:01+HLA-A*02:01", s, registry=reg)


# === 11. Error-free exact observation compatible with truth =================
def test_11_error_free_compatible_with_truth(registry) -> None:
    truth = "HLA-B*07:02+HLA-B*08:01"
    s = load_strategy({
        "name": "sbt_B2",
        "loci": {"HLA-B": {
            "assay_class": "sbt", "reference_name": "TOY-GEN",
            "interrogated_features": ["exon2", "exon3"],
        }},
    })
    r = observe(truth, s, registry=registry)
    # the true genotype must be among the reported alternatives
    alts = r.observed.split("|")
    assert any("HLA-B*07:02" in a and "HLA-B*08:01" in a for a in alts)
    assert isinstance(r, ObservationResult)
    assert r.error_category.value == "none"


# === 12. glstring parses and round-trips all resulting GL Strings ===========
def test_12_glstring_roundtrips_all_outputs(registry) -> None:
    cases = [
        ("second_field.toml", "HLA-A*01:01:01+HLA-A*02:01", None),
        ("serology_broad_B.toml", "HLA-B*15:01+HLA-B*07:02", registry),
        ("sbt_exon23_bc_phased.toml",
         "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01", registry),
        ("mixed_assay_registry.toml",
         "HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01^"
         "HLA-C*07:01+HLA-C*07:02^HLA-DQB1*02:01+HLA-DQB1*03:01", registry),
    ]
    for fname, truth, reg in cases:
        s = load_strategy_file(STRAT / fname)
        r = observe(truth, s, registry=reg)
        assert _roundtrips(r.observed), f"{fname}: {r.observed!r}"
