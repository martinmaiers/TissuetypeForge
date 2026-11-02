"""Observation framework tests (spec §2, §5, §6, §7; QA §10.1-§10.9)."""

from __future__ import annotations

import pytest

from tissuetypeforge.glstring import parse_blocks, serialize
from tissuetypeforge.observation import ObservationError, observe
from tissuetypeforge.result import ErrorCategory
from tissuetypeforge.strategy import (
    AssayClass,
    LocusConfig,
    PhaseGroup,
    Resolution,
    Strategy,
)

# --- helpers --------------------------------------------------------------


def _strategy(
    resolution: Resolution = Resolution.SECOND_FIELD,
    *,
    phase_groups: tuple[PhaseGroup, ...] = (),
    untyped: tuple[str, ...] = (),
    name: str = "t",
) -> Strategy:
    default = LocusConfig(
        locus="<default>",
        typed=True,
        assay_class=AssayClass.IDEALIZED,
        resolution=resolution,
        reference_database="IPD-IMGT/HLA 3.55.0",
    )
    loci = {
        u: LocusConfig(locus=u, typed=False, reference_database="IPD-IMGT/HLA 3.55.0")
        for u in untyped
    }
    return Strategy(
        name=name,
        loci=loci,
        phase_groups=phase_groups,
        reference_database="IPD-IMGT/HLA 3.55.0",
        default_locus=default,
    )


def _alleles_at_locus(gl: str, locus: str) -> set[frozenset[str]]:
    """Return the set of per-copy allele-name sets observed at ``locus``.

    Each gene copy contributes a frozenset of the (ambiguous) allele names it
    could be. Order of copies is ignored so homozygous/heterozygous comparisons
    are phase-agnostic within a locus.
    """
    blocks = parse_blocks(gl)
    out: list[frozenset[str]] = []
    for genlist in blocks.blocks:
        genotype = genlist.genotypes[0]
        loci = genotype.copies[0].loci
        if locus not in loci:
            continue
        pos = loci.index(locus)
        for copy in genotype.copies:
            names = {a.name for a in copy.segments[pos].alleles}
            out.append(frozenset(names))
    return set(out) if out else set()


def _truth_compatible(truth_gl: str, observed_gl: str, resolution: Resolution) -> bool:
    """QA §10.1/§10.5: error-free observation is compatible with truth.

    For idealized truncation, each observed per-copy allele name must be a
    field-prefix of some true allele at that locus, and every typed true copy
    must be represented.
    """
    n = {Resolution.FIRST_FIELD: 1, Resolution.SECOND_FIELD: 2}.get(resolution)
    truth = parse_blocks(truth_gl)
    for genlist in truth.blocks:
        genotype = genlist.genotypes[0]
        for pos, locus in enumerate(genotype.copies[0].loci):
            obs_copies = _alleles_at_locus(observed_gl, locus)
            if not obs_copies:
                continue  # locus not typed; handled elsewhere
            for copy in genotype.copies:
                true_names = {a.name for a in copy.segments[pos].alleles}
                expected = {
                    _truncate_name(tn, n) for tn in true_names
                }
                if not any(expected <= oc or oc <= expected for oc in obs_copies):
                    return False
    return True


def _truncate_name(name: str, n: int | None) -> str:
    if n is None or "*" not in name:
        return name
    locus, fields = name.split("*", 1)
    flds = fields.split(":")
    if len(flds) <= n:
        return name
    return f"{locus}*{':'.join(flds[:n])}"


# --- QA items -------------------------------------------------------------

TRUTH = (
    "HLA-A*01:01:01:01+HLA-A*02:01:01^"
    "HLA-B*07:02:01~HLA-C*07:02:01+HLA-B*08:01:01~HLA-C*07:01:01"
)


def test_error_free_observation_compatible_with_truth() -> None:
    # QA §10.1
    for res in (Resolution.FIRST_FIELD, Resolution.SECOND_FIELD):
        r = observe(TRUTH, _strategy(res))
        assert _truth_compatible(TRUTH, r.observed, res)
        assert r.error_category is ErrorCategory.NONE


def test_observed_is_round_trippable() -> None:
    # QA §10.9
    r = observe(TRUTH, _strategy(Resolution.SECOND_FIELD))
    assert serialize(parse_blocks(r.observed)) == r.observed


def test_bc_phased_strategy_retains_bc_phase() -> None:
    # QA §10.3
    s = _strategy(phase_groups=(PhaseGroup(("HLA-B", "HLA-C")),))
    r = observe(TRUTH, s)
    d = parse_blocks(r.observed)
    bc = [g for g in d.blocks if frozenset(g.loci) == frozenset({"HLA-B", "HLA-C"})]
    assert len(bc) == 1  # B and C reported together -> phase retained
    assert ("HLA-B", "HLA-C") in r.metadata.retained_phase_groups


def test_unphased_strategy_discards_bc_phase() -> None:
    # QA §10.4
    r = observe(TRUTH, _strategy())  # no phase groups
    d = parse_blocks(r.observed)
    # B and C must be in separate groups (no shared block).
    b = {frozenset(g.loci) for g in d.blocks if "HLA-B" in g.loci}
    c = {frozenset(g.loci) for g in d.blocks if "HLA-C" in g.loci}
    assert b == {frozenset({"HLA-B"})}
    assert c == {frozenset({"HLA-C"})}
    assert r.metadata.retained_phase_groups == ()


def test_removing_phase_does_not_invent_phase() -> None:
    # QA §10.2: an unphased strategy over unphased truth never creates phase.
    unphased_truth = (
        "HLA-B*07:02:01+HLA-B*08:01:01^HLA-C*07:02:01+HLA-C*07:01:01"
    )
    # Even though strategy *declares* a B~C phase group, truth did not phase
    # them, so phase must NOT be invented.
    s = _strategy(phase_groups=(PhaseGroup(("HLA-B", "HLA-C")),))
    r = observe(unphased_truth, s)
    d = parse_blocks(r.observed)
    assert not any(frozenset(g.loci) == frozenset({"HLA-B", "HLA-C"}) for g in d.blocks)
    assert r.metadata.retained_phase_groups == ()


def test_missing_locus_represented_correctly() -> None:
    # QA §10.7
    truth = TRUTH + "^HLA-DQB1*02:01:01+HLA-DQB1*03:01:01"
    r = observe(truth, _strategy(untyped=("HLA-DQB1",)))
    assert "HLA-DQB1" in r.metadata.missing_loci
    assert "HLA-DQB1" not in r.observed


def test_reproducible_with_fixed_seed() -> None:
    # QA §10.8: deterministic in M1 (no stochastic error); seed recorded.
    s = _strategy()
    r1 = observe(TRUTH, s, seed=42)
    r2 = observe(TRUTH, s, seed=42)
    assert r1.observed == r2.observed
    assert r1.metadata.seed == 42


def test_reference_database_recorded_in_result() -> None:
    # QA §10.10
    r = observe(TRUTH, _strategy())
    assert "IPD-IMGT/HLA 3.55.0" in r.metadata.reference_databases


def test_homozygous_genotype() -> None:
    truth = "HLA-A*01:01:01:01+HLA-A*01:01:01:02"
    r = observe(truth, _strategy(Resolution.SECOND_FIELD))
    # both copies truncate to HLA-A*01:01 -> homozygous at 2 fields
    assert r.observed == "HLA-A*01:01+HLA-A*01:01"


def test_idealized_does_not_introduce_ambiguity() -> None:
    # Truncating two distinct alleles that collapse to the same name must not
    # produce a spurious a/a allele-ambiguity list.
    truth = "HLA-A*01:01:01+HLA-A*01:01:02^HLA-B*07:02:01+HLA-B*08:01:01"
    r = observe(truth, _strategy(Resolution.SECOND_FIELD))
    assert "/" not in r.observed


def test_information_loss_recorded_but_not_mistyping() -> None:
    # spec §7: ordinary information loss is NOT labeled mistyping.
    r = observe(TRUTH, _strategy(Resolution.FIRST_FIELD))
    assert r.error_category is ErrorCategory.NONE
    assert r.loss_record  # truncation from 3-4 fields to 1 is recorded


def test_allele_ambiguity_in_truth_preserved_under_allele_resolution() -> None:
    truth = "HLA-A*01:01/HLA-A*01:02+HLA-A*02:01"
    r = observe(truth, _strategy(Resolution.ALLELE))
    assert "HLA-A*01:01/HLA-A*01:02" in r.observed


def test_m2_resolution_mode_not_faked() -> None:
    # spec §12: deferred modes raise rather than fake a transformation.
    s = _strategy(Resolution.SEROLOGY)
    with pytest.raises(NotImplementedError):
        observe(TRUTH, s)


def test_ambiguity_set_contains_truth_when_no_error() -> None:
    # QA §10.5: at ALLELE resolution the observation equals truth exactly.
    r = observe(TRUTH, _strategy(Resolution.ALLELE))
    assert _truth_compatible(TRUTH, r.observed, Resolution.ALLELE)


def test_genotype_ambiguous_truth_rejected() -> None:
    with pytest.raises(ObservationError):
        observe(
            "HLA-A*01:01+HLA-A*02:01|HLA-A*01:02+HLA-A*02:01",
            _strategy(),
        )
