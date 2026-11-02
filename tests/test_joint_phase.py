"""M2B tests — joint-phase molecule enumeration (§6.4.1, D8).

Each locus individually ambiguous; the joint (whole-molecule) observation yields
the product of molecules (4), never the per-locus cross-product (8). The
fabricated trans diplotype B*A~C*B + B*2~C*1 must be absent.
"""

from __future__ import annotations

from tissuetypeforge.glstring import serialize
from tissuetypeforge.observation.joint_phase import (
    enumerate_joint,
    joint_to_genotype_list,
)
from tissuetypeforge.reference.model import Base, ColumnModel, FeatureKey

EX2 = FeatureKey("exon", 2)

# B locus: B*A and B*Aq share interrogated column 'G'; B*2 is 'T'.
_B = ColumnModel(
    locus="HLA-B",
    database_version="TOY/B",
    feature_of_column=(EX2,),
    rows={
        "HLA-B*A": (Base.G,),
        "HLA-B*Aq": (Base.G,),   # indistinguishable from B*A over interrogated col
        "HLA-B*2": (Base.T,),
    },
)
# C locus: C*B and C*Bq share interrogated column 'A'; C*1 is 'G'.
_C = ColumnModel(
    locus="HLA-C",
    database_version="TOY/C",
    feature_of_column=(EX2,),
    rows={
        "HLA-C*1": (Base.G,),
        "HLA-C*B": (Base.A,),
        "HLA-C*Bq": (Base.A,),   # indistinguishable from C*B over interrogated col
    },
)

_MODELS = {"HLA-B": _B, "HLA-C": _C}
_LOCI = ("HLA-B", "HLA-C")
_FEATURES = {"HLA-B": frozenset({EX2}), "HLA-C": frozenset({EX2})}
# truth molecules: hap1 = B*A ~ C*1 ; hap2 = B*2 ~ C*B
_TRUTH = (("HLA-B*A", "HLA-C*1"), ("HLA-B*2", "HLA-C*B"))


def test_joint_phase_product_of_molecules_not_loci() -> None:
    cand = enumerate_joint(_MODELS, _LOCI, _TRUTH, _FEATURES)
    # Exactly 4 candidate diplotypes (2 B-class x 2 C-class molecules), NOT 8.
    assert len(cand.confirmed) == 4
    pairs = set(cand.confirmed)
    # truth is present
    assert tuple(sorted((("HLA-B*A", "HLA-C*1"), ("HLA-B*2", "HLA-C*B")))) in pairs
    # fabricated trans diplotype must be ABSENT
    fabricated = tuple(
        sorted((("HLA-B*A", "HLA-C*B"), ("HLA-B*2", "HLA-C*1")))
    )
    assert fabricated not in pairs


def test_joint_phase_serializes_phased_haplotypes() -> None:
    cand = enumerate_joint(_MODELS, _LOCI, _TRUTH, _FEATURES)
    gl = serialize(joint_to_genotype_list(cand))
    # every alternative must phase B with C via '~', and there are 4 alternatives
    assert gl.count("|") == 3  # 4 alternatives
    assert "~" in gl
    # the fabricated cis haplotype B*A~C*B must not appear
    assert "HLA-B*A~HLA-C*B" not in gl


def test_joint_phase_no_fabricated_haplotype_in_any_alternative() -> None:
    cand = enumerate_joint(_MODELS, _LOCI, _TRUTH, _FEATURES)
    # C*1 only ever pairs (cis) with a B*A-class allele; C*B/C*Bq only with B*2.
    for h1, h2 in cand.confirmed:
        for hap in (h1, h2):
            b, c = hap
            if c == "HLA-C*1":
                assert b in {"HLA-B*A", "HLA-B*Aq"}
            else:  # C*B or C*Bq
                assert b == "HLA-B*2"
