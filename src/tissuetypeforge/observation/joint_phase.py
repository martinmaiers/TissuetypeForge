"""M2B — joint-phase (whole-molecule) enumeration for a declared group (D8).

When a strategy declares a ``joint_observation_group`` (member loci + the shared
interrogated region) backed by a phase-resolving assay, candidate genotypes are
enumerated as **haplotype (molecule) pairings over the concatenated joint
signature** — the product of observed molecules, NOT the product of per-locus
ambiguities (main §6.4.1). This prevents inventing cross-locus haplotypes.

Scope note: this reference implementation handles a two-molecule (diploid) group
whose member loci each contribute one interrogated column-set, modeled as a
per-locus :class:`ColumnModel`. Truth is the two phased haplotypes.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement, product

from tissuetypeforge.glstring import (
    Allele,
    AlleleList,
    Genotype,
    GenotypeList,
    PhasedGenotype,
)
from tissuetypeforge.observation.policy import (
    DEFAULT_LIMITS,
    IndeterminateReferenceError,
    ReferencePolicy,
    ResourceLimits,
)
from tissuetypeforge.reference.model import Base, ColumnModel
from tissuetypeforge.reference.signature import (
    ColumnComparison,
    interrogated_columns,
    signature,
)

# A haplotype over the ordered group loci: one allele name per locus.
Haplotype = tuple[str, ...]


@dataclass(frozen=True)
class JointCandidateSet:
    loci: tuple[str, ...]
    #: Each candidate diplotype is an unordered pair of haplotypes.
    confirmed: tuple[tuple[Haplotype, Haplotype], ...]
    indeterminate: tuple[tuple[Haplotype, Haplotype], ...]
    database_versions: tuple[str, ...]
    indeterminate_alleles: tuple[str, ...]


def _hap_signature(
    models: dict[str, ColumnModel],
    cols: dict[str, tuple[int, ...]],
    hap: Haplotype,
    loci: tuple[str, ...],
) -> tuple[Base, ...] | None:
    """Concatenated signature of a haplotype across the group's loci."""
    out: list[Base] = []
    for locus, allele in zip(loci, hap, strict=True):
        s = signature(models[locus], allele, cols[locus])
        if s is None:
            return None
        out.extend(s)
    return tuple(out)


def _compare_molecule_multiset(
    observed: tuple[tuple[Base, ...], tuple[Base, ...]],
    cand: tuple[tuple[Base, ...], tuple[Base, ...]],
) -> ColumnComparison:
    """Compare the observed pair of molecule signatures to a candidate pair.

    Molecules are unordered, so try both assignments; a candidate matches if
    some assignment is not INCOMPATIBLE. INDETERMINATE if a gap/missing base is
    involved in the best assignment.
    """
    best = ColumnComparison.INCOMPATIBLE
    for cand_ordered in (cand, (cand[1], cand[0])):
        res = ColumnComparison.COMPATIBLE
        for obs_sig, cand_sig in zip(observed, cand_ordered, strict=True):
            for ob, cb in zip(obs_sig, cand_sig, strict=True):
                if ob is Base.INDETERMINATE or cb is Base.INDETERMINATE:
                    res = ColumnComparison.INDETERMINATE
                elif ob != cb:
                    res = ColumnComparison.INCOMPATIBLE
                    break
            if res is ColumnComparison.INCOMPATIBLE:
                break
        # keep the best (COMPATIBLE > INDETERMINATE > INCOMPATIBLE)
        if res is ColumnComparison.COMPATIBLE:
            return ColumnComparison.COMPATIBLE
        if res is ColumnComparison.INDETERMINATE:
            best = ColumnComparison.INDETERMINATE
    return best


def enumerate_joint(
    models: dict[str, ColumnModel],
    loci: tuple[str, ...],
    true_haplotypes: tuple[Haplotype, Haplotype],
    interrogated_features: dict[str, frozenset],
    *,
    reference_policy: ReferencePolicy = ReferencePolicy.STRICT,
    limits: ResourceLimits = DEFAULT_LIMITS,
) -> JointCandidateSet:
    """Enumerate candidate diplotypes as molecule pairings (main §6.4.1, D8)."""
    cols = {
        locus: interrogated_columns(models[locus], interrogated_features[locus])
        for locus in loci
    }
    # observed molecule signatures (truth's two haplotypes)
    obs = tuple(_hap_signature(models, cols, h, loci) for h in true_haplotypes)
    if any(s is None for s in obs):
        raise IndeterminateReferenceError(
            "true haplotype(s) absent from reference over the joint region",
            indeterminate_alleles=(),
        )
    observed: tuple[tuple[Base, ...], tuple[Base, ...]] = (obs[0], obs[1])  # type: ignore[assignment]

    # candidate haplotype universe = product of per-locus alleles with a signature
    per_locus_alleles: list[list[str]] = []
    for locus in loci:
        usable = [a for a in models[locus].alleles()
                  if signature(models[locus], a, cols[locus]) is not None]
        per_locus_alleles.append(usable)
    haps: list[Haplotype] = [tuple(h) for h in product(*per_locus_alleles)]
    limits.check_universe(len(haps))

    hap_sig: dict[Haplotype, tuple[Base, ...]] = {}
    for h in haps:
        s = _hap_signature(models, cols, h, loci)
        if s is not None:
            hap_sig[h] = s

    confirmed: list[tuple[Haplotype, Haplotype]] = []
    indeterminate: list[tuple[Haplotype, Haplotype]] = []
    indet_alleles: set[str] = set()

    count = 0
    for h1, h2 in combinations_with_replacement(sorted(hap_sig), 2):
        count += 1
        cmp = _compare_molecule_multiset(observed, (hap_sig[h1], hap_sig[h2]))
        if cmp is ColumnComparison.INCOMPATIBLE:
            continue
        pair = tuple(sorted((h1, h2)))
        if cmp is ColumnComparison.COMPATIBLE:
            confirmed.append(pair)  # type: ignore[arg-type]
        else:
            indeterminate.append(pair)  # type: ignore[arg-type]
            indet_alleles.update(h1)
            indet_alleles.update(h2)
    limits.check_joint_molecules(count)

    if indeterminate and reference_policy is ReferencePolicy.STRICT:
        raise IndeterminateReferenceError(
            f"strict reference_policy: incomplete reference leaves "
            f"{len(indeterminate)} joint candidate(s) undecidable over {loci}",
            indeterminate_alleles=tuple(sorted(indet_alleles)),
        )

    return JointCandidateSet(
        loci=loci,
        confirmed=tuple(sorted(set(confirmed))),
        indeterminate=tuple(sorted(set(indeterminate))),
        database_versions=tuple(sorted({models[x].database_version for x in loci})),
        indeterminate_alleles=tuple(sorted(indet_alleles)),
    )


def _haplotype_to_phased(hap: Haplotype) -> PhasedGenotype:
    return PhasedGenotype(tuple(AlleleList((Allele(a),)) for a in hap))


def joint_to_genotype_list(cand: JointCandidateSet, *,
                           limits: ResourceLimits = DEFAULT_LIMITS) -> GenotypeList:
    """Serialize joint candidates to one phased ``GenotypeList`` over the group.

    Each alternative is a ``Genotype`` of two phased haplotype copies (``~``
    within a copy, ``+`` between copies). Alternatives joined by ``|`` are the
    exact molecule-pairing candidates (no cross-locus invention).
    """
    pairs = tuple(sorted(set(cand.confirmed) | set(cand.indeterminate)))
    if not pairs:
        raise IndeterminateReferenceError(
            f"no joint candidate genotype over {cand.loci}",
            indeterminate_alleles=cand.indeterminate_alleles,
        )
    limits.check_alternatives(len(pairs))
    genotypes = tuple(
        Genotype((_haplotype_to_phased(h1), _haplotype_to_phased(h2)))
        for h1, h2 in pairs
    )
    return GenotypeList(genotypes)
