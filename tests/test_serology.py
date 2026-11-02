"""M2C tests — serology backend (verified broad/split worked examples).

Toy artifact (companion §6):
  B*07:02 -> {B7}     B*08:01 -> {B8}
  B*15:01 -> {B62}    B*15:02 -> {B75}   B*15:03 -> {B72}  B*15:99 -> {B15}
  lattice: B62,B75,B72,B15 broad_of -> B15 ; B7->B7 ; B8->B8
"""

from __future__ import annotations

from tissuetypeforge.glstring import serialize
from tissuetypeforge.serology import (
    Assignment,
    BroadSplitLattice,
    Confidence,
    Epitope,
    ResolutionLevel,
    SerologyModel,
    Specificity,
    observe_serology,
)

B7, B8 = Specificity("B7"), Specificity("B8")
B62, B75, B72, B15 = (
    Specificity("B62"),
    Specificity("B75"),
    Specificity("B72"),
    Specificity("B15"),
)
Bw4, Bw6 = Epitope("Bw4"), Epitope("Bw6")

_LATTICE = BroadSplitLattice(
    broad_of={B62: B15, B75: B15, B72: B15, B15: B15, B7: B7, B8: B8}
)


def _asg(specs, conf=Confidence.ASSIGNED, epitopes=frozenset()):
    return Assignment(frozenset(specs), conf, epitopes)


_MODEL = SerologyModel(
    locus="HLA-B",
    database_version="TOY-SER/1.0",
    assignment={
        "HLA-B*07:02": _asg({B7}, epitopes={Bw6}),
        "HLA-B*08:01": _asg({B8}, epitopes={Bw6}),
        "HLA-B*15:01": _asg({B62}, epitopes={Bw6}),
        "HLA-B*15:02": _asg({B75}, epitopes={Bw4}),
        "HLA-B*15:03": _asg({B72}, epitopes={Bw6}),
        "HLA-B*15:99": _asg({B15}, epitopes={Bw6}),
    },
    lattice=_LATTICE,
)


def test_split_level_unambiguous() -> None:
    # §6.1: B*15:01 + B*07:02 at split -> {B62,B7}, single compatible genotype.
    obs = observe_serology(_MODEL, ("HLA-B*15:01", "HLA-B*07:02"), ResolutionLevel.SPLIT)
    assert obs.report.specificities == frozenset({B62, B7})
    assert obs.report.resolution_level is ResolutionLevel.SPLIT
    assert obs.compatible_pairs == (("HLA-B*07:02", "HLA-B*15:01"),)
    assert serialize(obs.compatible_genotypes) == "HLA-B*07:02+HLA-B*15:01"


def test_broad_only_four_way_ambiguity() -> None:
    # §6.2: broad -> {B15,B7}; four compatible DNA pairs; factors to '/'.
    obs = observe_serology(_MODEL, ("HLA-B*15:01", "HLA-B*07:02"), ResolutionLevel.BROAD)
    assert obs.report.specificities == frozenset({B15, B7})
    assert set(obs.compatible_pairs) == {
        ("HLA-B*07:02", "HLA-B*15:01"),
        ("HLA-B*07:02", "HLA-B*15:02"),
        ("HLA-B*07:02", "HLA-B*15:03"),
        ("HLA-B*07:02", "HLA-B*15:99"),
    }
    gl = serialize(obs.compatible_genotypes)
    assert (
        gl
        == "HLA-B*07:02+HLA-B*15:01/HLA-B*15:02/HLA-B*15:03/HLA-B*15:99"
    )


def test_broad_only_never_emits_split_tokens() -> None:
    obs = observe_serology(_MODEL, ("HLA-B*15:01", "HLA-B*07:02"), ResolutionLevel.BROAD)
    assert B62 not in obs.report.specificities
    assert B75 not in obs.report.specificities
    assert obs.report.specificities == frozenset({B15, B7})


def test_antigen_tokens_never_in_compatible_gl_string() -> None:
    obs = observe_serology(_MODEL, ("HLA-B*15:01", "HLA-B*07:02"), ResolutionLevel.BROAD)
    gl = serialize(obs.compatible_genotypes)
    for tok in ("B7", "B15", "B62", "Bw4", "Bw6"):
        # the only allowed occurrences are inside DNA names like HLA-B*15:01
        assert f"+{tok}" not in gl and f"/{tok}" not in gl and not gl.startswith(tok)


def test_is_not_name_truncation() -> None:
    # B*15:01 serologic specificity is B62 (via artifact), never 'B15' from name.
    obs = observe_serology(_MODEL, ("HLA-B*15:01", "HLA-B*07:02"), ResolutionLevel.SPLIT)
    assert B62 in obs.report.specificities
    assert B15 not in obs.report.specificities


def test_bw4_bw6_reported_on_separate_axis() -> None:
    obs = observe_serology(
        _MODEL, ("HLA-B*15:02", "HLA-B*07:02"), ResolutionLevel.SPLIT,
        observed_epitopes=frozenset({Bw4, Bw6}),
    )
    # B*15:02 carries Bw4, B*07:02 carries Bw6
    assert obs.report.public_epitopes == frozenset({Bw4, Bw6})

def test_no_independent_cross_product() -> None:
    # §6.6: cross-reactive allele c -> {B7,B15}; observed broad {B7,B15}.
    # A naive a/b/c + a/b/c union cross-product would assert excluded a+a ({B7})
    # and b+b ({B15}); the engine must enumerate the exact pair set instead.
    model = SerologyModel(
        locus="HLA-B",
        database_version="TOY-SER/66",
        assignment={
            "HLA-B*07:02": _asg({B7}),          # a -> B7
            "HLA-B*08:01": _asg({B15}),         # b -> B15
            "HLA-B*15:16": _asg({B7, B15}),     # c -> cross-reactive
        },
        lattice=BroadSplitLattice(broad_of={B7: B7, B15: B15}),
    )
    obs = observe_serology(model, ("HLA-B*07:02", "HLA-B*08:01"), ResolutionLevel.BROAD)
    pairs = set(obs.compatible_pairs)
    assert pairs == {
        ("HLA-B*07:02", "HLA-B*08:01"),
        ("HLA-B*07:02", "HLA-B*15:16"),
        ("HLA-B*08:01", "HLA-B*15:16"),
        ("HLA-B*15:16", "HLA-B*15:16"),
    }
    assert ("HLA-B*07:02", "HLA-B*07:02") not in pairs  # excluded a+a
    assert ("HLA-B*08:01", "HLA-B*08:01") not in pairs  # excluded b+b



def test_null_allele_is_blank_from_not_expressed() -> None:
    model = SerologyModel(
        locus="HLA-B",
        database_version="TOY-SER/null",
        assignment={
            "HLA-B*07:02": _asg({B7}),
            "HLA-B*15:01:01:02N": _asg(frozenset(), conf=Confidence.NOT_EXPRESSED),
        },
        lattice=_LATTICE,
    )
    obs = observe_serology(
        model, ("HLA-B*07:02", "HLA-B*15:01:01:02N"), ResolutionLevel.SPLIT
    )
    # null contributes blank -> observed is just {B7}
    assert obs.report.specificities == frozenset({B7})
    assert obs.report.blank_present is True
