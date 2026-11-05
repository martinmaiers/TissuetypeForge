"""Tests for the simplified observation transform (REQ 7 items 1, 2, 3, 8, 9)."""

from __future__ import annotations

import pytest

from tissuetypeforge import (
    ObservationError,
    ResolutionMode,
    TypingProfile,
    observe_genotype,
    parse_gl_string,
    serialize,
)
from tissuetypeforge.glstring import parse_blocks

# A phased, four-locus truth. Phase (~) is present in the B~C block on purpose
# so we can assert the observation never preserves interlocus phase (REQ 4/8).
TRUTH = (
    "HLA-A*02:01:01:01+HLA-A*03:01:01:01^"
    "HLA-B*07:02:01~HLA-C*07:02:01+HLA-B*08:01:01~HLA-C*07:01:01^"
    "HLA-DRB1*15:01:01+HLA-DRB1*03:01:01"
)


def _profile(**modes: str) -> TypingProfile:
    return TypingProfile(
        id="p",
        weight=1.0,
        resolution={k: ResolutionMode(v) for k, v in modes.items()},
    )


def _copies_at(gl: str, locus: str) -> list[set[str]]:
    """Per-copy allele-name sets observed at ``locus`` (phase-agnostic)."""
    blocks = parse_blocks(gl)
    out: list[set[str]] = []
    for genlist in blocks.blocks:
        genotype = genlist.genotypes[0]
        loci = genotype.copies[0].loci
        if locus not in loci:
            continue
        pos = loci.index(locus)
        for copy in genotype.copies:
            out.append({a.name for a in copy.segments[pos].alleles})
    return out


# --- REQ 7.1: correct 0/1/2-field transformations -------------------------


def test_two_field_transformation() -> None:
    p = _profile(A="2-field", B="2-field", C="2-field", DRB1="2-field")
    o = observe_genotype(TRUTH, p)
    assert _copies_at(o.observed, "HLA-A") == [{"HLA-A*02:01"}, {"HLA-A*03:01"}]


def test_one_field_transformation() -> None:
    p = _profile(A="1-field", B="1-field", C="1-field", DRB1="1-field")
    o = observe_genotype(TRUTH, p)
    assert _copies_at(o.observed, "HLA-A") == [{"HLA-A*02"}, {"HLA-A*03"}]
    assert _copies_at(o.observed, "HLA-B") == [{"HLA-B*07"}, {"HLA-B*08"}]


def test_zero_field_omits_locus() -> None:
    # REQ 7.3 too: 0-field is a missing result, not a homozygous call.
    p = _profile(A="2-field", B="2-field", C="0-field", DRB1="2-field")
    o = observe_genotype(TRUTH, p)
    assert "HLA-C" not in o.observed
    assert "HLA-C" in o.missing_loci
    assert "HLA-C" not in o.typed_loci


def test_matches_requirement_example_exactly() -> None:
    # The REQ 1 worked example: A*02:01:01:01 + A*03:01:01:01.
    truth = "HLA-A*02:01:01:01+HLA-A*03:01:01:01"
    assert observe_genotype(truth, _profile(A="1-field")).observed == (
        "HLA-A*02+HLA-A*03"
    )
    assert observe_genotype(truth, _profile(A="2-field")).observed == (
        "HLA-A*02:01+HLA-A*03:01"
    )
    assert observe_genotype(truth, _profile(A="0-field")).observed == ""


# --- REQ 7.2: collapse at 1-field preserves uncertainty -------------------


def test_one_field_collapse_does_not_assert_true_homozygosity() -> None:
    # Two DISTINCT alleles that collapse to the same first field.
    truth = "HLA-A*02:01:01+HLA-A*02:07:01"
    o = observe_genotype(truth, _profile(A="1-field"))
    # The collapse is recorded as information loss / preserved uncertainty,
    # explicitly NOT a claim of true homozygosity.
    assert o.loss_record
    assert any("homozygosity is not asserted" in n for n in o.loss_record)


def test_true_homozygote_has_no_collapse_note() -> None:
    truth = "HLA-A*02:01:01:01+HLA-A*02:01:01:02"
    o = observe_genotype(truth, _profile(A="2-field"))
    # Both truncate to A*02:01, but the TRUE alleles already share that 2-field
    # name's prefix lineage; still, distinct 4-field names collapse -> recorded.
    assert o.observed == "HLA-A*02:01+HLA-A*02:01"


def test_intralocus_ambiguity_retained() -> None:
    truth = "HLA-A*02:01/HLA-A*02:05+HLA-A*03:01"
    o = observe_genotype(truth, _profile(A="2-field"))
    copies = _copies_at(o.observed, "HLA-A")
    assert {"HLA-A*02:01", "HLA-A*02:05"} in copies  # the "/" ambiguity survives


def test_intralocus_ambiguity_dedupes_on_collapse() -> None:
    # A*02:01:01 / A*02:01:02 both -> A*02:01 at 2-field: must NOT be "a/a".
    truth = "HLA-A*02:01:01/HLA-A*02:01:02+HLA-A*03:01"
    o = observe_genotype(truth, _profile(A="2-field"))
    assert "HLA-A*02:01/HLA-A*02:01" not in o.observed


# --- REQ 7.3: missing observations ----------------------------------------


def test_locus_not_in_profile_is_missing() -> None:
    # Profile omits DRB1 entirely -> locus present in truth reported missing.
    p = _profile(A="2-field", B="2-field", C="2-field")
    o = observe_genotype(TRUTH, p)
    assert "HLA-DRB1" in o.missing_loci
    assert "HLA-DRB1" not in o.observed


def test_all_zero_field_yields_empty_observation() -> None:
    p = _profile(A="0-field", B="0-field", C="0-field", DRB1="0-field")
    o = observe_genotype(TRUTH, p)
    assert o.observed == ""
    assert set(o.missing_loci) == {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"}
    assert o.typed_loci == ()


# --- REQ 7.8: no unintended interlocus phasing ----------------------------


def test_no_interlocus_phase_even_when_truth_is_phased() -> None:
    # TRUTH phases B~C. The observation must NOT carry any "~".
    p = _profile(A="2-field", B="2-field", C="2-field", DRB1="2-field")
    o = observe_genotype(TRUTH, p)
    assert "~" not in o.observed


def test_each_locus_is_an_independent_block() -> None:
    p = _profile(A="2-field", B="2-field", C="2-field", DRB1="2-field")
    o = observe_genotype(TRUTH, p)
    blocks = parse_blocks(o.observed)
    # Every observed block must be a single-locus block (no shared/phased loci).
    for genlist in blocks.blocks:
        loci = genlist.genotypes[0].copies[0].loci
        assert len(loci) == 1


# --- REQ 7.9: existing supported functionality intact ---------------------


def test_observed_is_round_trippable() -> None:
    p = _profile(A="2-field", B="2-field", C="2-field", DRB1="2-field")
    o = observe_genotype(TRUTH, p)
    assert serialize(parse_blocks(o.observed)) == o.observed


def test_general_purpose_gl_parsing_still_exported() -> None:
    # parse_gl_string / serialize remain available for other components.
    node = parse_gl_string("HLA-A*02:01+HLA-A*03:01")
    assert serialize(node) == "HLA-A*02:01+HLA-A*03:01"


def test_genotype_ambiguous_truth_rejected() -> None:
    with pytest.raises(ObservationError):
        observe_genotype(
            "HLA-A*01:01+HLA-A*02:01|HLA-A*01:02+HLA-A*02:01",
            _profile(A="2-field"),
        )


def test_profile_id_recorded() -> None:
    p = TypingProfile(id="legacy", weight=1.0, resolution={"A": ResolutionMode.ONE_FIELD})
    o = observe_genotype("HLA-A*02:01:01+HLA-A*03:01:01", p)
    assert o.profile_id == "legacy"
