"""M2C — serology observation backend (serology companion design).

Produces two separate outputs (companion §7):

1. a :class:`SerologicAntigenReport` (antigen/epitope tokens; NOT a GL String);
2. a compatible-DNA genotype set as a GL ``GenotypeList`` of real DNA alleles.

The genotype-ambiguity algorithm is the same pair-enumeration skeleton as the
sequence backend: enumerate unordered DNA pairs whose antigen observation
reproduces truth's, then reduce with exact ``|`` / lossless ``/`` (companion §6).
A broad-only assay never emits split-level tokens (companion §6.2).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement

from tissuetypeforge.glstring import GenotypeList
from tissuetypeforge.observation.policy import DEFAULT_LIMITS, ResourceLimits
from tissuetypeforge.observation.reduce import pairs_to_genotype_list
from tissuetypeforge.serology.model import (
    Confidence,
    Epitope,
    ResolutionLevel,
    SerologicAntigenReport,
    SerologyModel,
    Specificity,
)

Pair = tuple[str, str]


@dataclass(frozen=True)
class SerologyObservation:
    """Combined result: the antigen report plus the compatible-DNA genotype set."""

    report: SerologicAntigenReport
    compatible_genotypes: GenotypeList
    compatible_pairs: tuple[Pair, ...]


def _observed_specificities(
    model: SerologyModel, a1: str, a2: str, level: ResolutionLevel
) -> frozenset[Specificity]:
    s1 = model.projected_specificities(a1, level) or frozenset()
    s2 = model.projected_specificities(a2, level) or frozenset()
    return s1 | s2


def _observed_epitopes(model: SerologyModel, a1: str, a2: str) -> frozenset[Epitope]:
    out: set[Epitope] = set()
    for a in (a1, a2):
        asg = model.get(a)
        if asg is not None:
            out |= asg.epitopes
    return out


def observe_serology(
    model: SerologyModel,
    true_pair: Pair,
    level: ResolutionLevel,
    *,
    observed_epitopes: frozenset[Epitope] = frozenset(),
    limits: ResourceLimits = DEFAULT_LIMITS,
) -> SerologyObservation:
    """Observe a serologic typing of ``true_pair`` at ``level`` (companion §4-§7).

    ``observed_epitopes`` is the explicitly declared set of public epitopes the
    assay interrogates (Decision 1). Only epitopes in this set may participate in
    the simulated observation: they are reported *and* constrain the compatible
    genotype set. Epitopes carried by alleles in the reference model but NOT in
    ``observed_epitopes`` are never observed and never narrow candidates — the
    reference-model assignment and the simulated observation are kept distinct.
    An empty set (default) means the assay does not interrogate any epitope.
    """
    a1, a2 = true_pair
    observed_spec = _observed_specificities(model, a1, a2, level)
    # Restrict epitope observation to those the assay declares it interrogates.
    observed_epi = _observed_epitopes(model, a1, a2) & observed_epitopes

    # confidence flags for the report
    flags = {"possible": False, "unknown": False, "blank": False}
    for a in (a1, a2):
        asg = model.get(a)
        if asg is None:
            continue
        if asg.confidence is Confidence.POSSIBLE:
            flags["possible"] = True
        elif asg.confidence is Confidence.UNKNOWN:
            flags["unknown"] = True
        elif asg.confidence is Confidence.NOT_EXPRESSED:
            flags["blank"] = True

    report = SerologicAntigenReport(
        locus=model.locus,
        specificities=observed_spec,  # at `level`; broad-only never has splits
        resolution_level=level,
        public_epitopes=observed_epi,  # only declared+observed epitopes
        database_version=model.database_version,
        uses_possible_assignment=flags["possible"],
        uses_unknown_assignment=flags["unknown"],
        blank_present=flags["blank"],
    )

    # Enumerate compatible DNA pairs: union of projected specificities (and, only
    # for declared epitopes, their restricted epitope set) equals the observation.
    # Pair enumeration, never independent per-allele.
    universe = model.alleles()
    limits.check_universe(len(universe))
    compatible: list[Pair] = []
    count = 0
    for x, y in combinations_with_replacement(universe, 2):
        count += 1
        if _observed_specificities(model, x, y, level) != observed_spec:
            continue
        if observed_epitopes and (
            _observed_epitopes(model, x, y) & observed_epitopes
        ) != observed_epi:
            continue
        compatible.append((x, y) if x <= y else (y, x))
    limits.check_candidate_pairs(count)
    compatible = sorted(set(compatible))

    if not compatible:
        raise ValueError("no compatible DNA genotype for the observed antigens")
    limits.check_alternatives(len(compatible))
    genotypes = pairs_to_genotype_list(compatible)
    return SerologyObservation(
        report=report,
        compatible_genotypes=genotypes,
        compatible_pairs=tuple(compatible),
    )
