"""M2B — diploid genotype enumeration for a single sequence-typed locus.

The scientific core of Milestone 2 (main design §3.2-§3.3):

1. Build the candidate universe ``U`` (all alleles in the pinned column model;
   incomplete-reference alleles stay in ``U`` — §4.2, never silently dropped).
2. Compute the true diploid *combined observation* ``O``: per interrogated
   column, the SET of bases across the two true alleles (§2.3; a set already
   encodes zygosity for a clean diploid).
3. Enumerate unordered candidate allele pairs whose combined observation is not
   positively INCOMPATIBLE with ``O`` (three-valued; §4.1).
4. Apply ``reference_policy`` to INDETERMINATE candidates (§4.0): ``strict``
   raises, ``permissive`` keeps-and-flags.
5. Reduce to the tightest *faithful* GL form — exact ``|`` by default, lossless
   ``/`` only when the candidate set equals the full Cartesian product (§3.3).

Independent per-allele matching is deliberately NOT used (that is the bug §7.2
demonstrates); enumeration is over pairs.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement

from tissuetypeforge.glstring import GenotypeList
from tissuetypeforge.observation.policy import (
    DEFAULT_LIMITS,
    IndeterminateReferenceError,
    MissingReferenceError,
    ReferencePolicy,
    RepresentationError,
    RepresentationPolicy,
    ResourceLimits,
)
from tissuetypeforge.observation.reduce import pairs_to_genotype_list
from tissuetypeforge.reference.model import Base, ColumnModel
from tissuetypeforge.reference.signature import (
    ColumnComparison,
    interrogated_columns,
    signature,
)

Pair = tuple[str, str]  # unordered, stored sorted


@dataclass(frozen=True)
class LocusCandidateSet:
    """Result of single-locus enumeration before GL serialization."""

    locus: str
    confirmed: tuple[Pair, ...]
    indeterminate: tuple[Pair, ...]
    database_version: str
    indeterminate_alleles: tuple[str, ...]

    @property
    def all_pairs(self) -> tuple[Pair, ...]:
        return tuple(sorted(set(self.confirmed) | set(self.indeterminate)))


def combined_observation(
    model: ColumnModel, a1: str, a2: str, columns: tuple[int, ...]
) -> tuple[frozenset[Base], ...]:
    """True diploid observation: per column, the set of the two alleles' bases.

    Raises :class:`MissingReferenceError` if either true allele is absent from
    the reference over the interrogated region (we cannot observe what we have
    no reference for — the true allele's own signature must exist).
    """
    s1 = signature(model, a1, columns)
    s2 = signature(model, a2, columns)
    if s1 is None or s2 is None:
        missing = [a for a, s in ((a1, s1), (a2, s2)) if s is None]
        raise MissingReferenceError(
            f"true allele(s) {missing} absent from reference {model.database_version}"
        )
    return tuple(frozenset({b1, b2}) for b1, b2 in zip(s1, s2, strict=True))


def _compare_pair(
    observed: tuple[frozenset[Base], ...],
    cand_sig_1: tuple[Base, ...],
    cand_sig_2: tuple[Base, ...],
) -> ColumnComparison:
    """Three-valued comparison of a candidate pair's observation to ``observed``.

    The candidate pair's own per-column observation is the set
    ``{base(x,c), base(y,c)}``; it is compared to the observed set. Aggregate:
    any INCOMPATIBLE column -> INCOMPATIBLE; else any INDETERMINATE -> INDETERMINATE;
    else COMPATIBLE.
    """
    result = ColumnComparison.COMPATIBLE
    for obs_col, b1, b2 in zip(observed, cand_sig_1, cand_sig_2, strict=True):
        cand_col = frozenset({b1, b2})
        col = _compare_column_sets(obs_col, cand_col)
        if col is ColumnComparison.INCOMPATIBLE:
            return ColumnComparison.INCOMPATIBLE
        if col is ColumnComparison.INDETERMINATE:
            result = ColumnComparison.INDETERMINATE
    return result


def _compare_column_sets(
    observed: frozenset[Base], candidate: frozenset[Base]
) -> ColumnComparison:
    """Compare observed base-set to a candidate pair's base-set at one column."""
    if Base.INDETERMINATE in candidate or Base.INDETERMINATE in observed:
        return ColumnComparison.INDETERMINATE
    if observed == candidate:
        return ColumnComparison.COMPATIBLE
    return ColumnComparison.INCOMPATIBLE


def enumerate_locus(
    model: ColumnModel,
    true_pair: Pair,
    interrogated_features: frozenset,
    *,
    reference_policy: ReferencePolicy = ReferencePolicy.STRICT,
    limits: ResourceLimits = DEFAULT_LIMITS,
) -> LocusCandidateSet:
    """Enumerate candidate diploid genotype pairs for one locus (main §3.2)."""
    columns = interrogated_columns(model, interrogated_features)
    a1, a2 = true_pair
    observed = combined_observation(model, a1, a2, columns)

    universe = model.alleles()
    limits.check_universe(len(universe))

    # Precompute candidate signatures once (bucketing is an optimization; the
    # reference implementation enumerates pairs directly for clarity).
    sigs: dict[str, tuple[Base, ...]] = {}
    for u in universe:
        s = signature(model, u, columns)
        if s is not None:
            sigs[u] = s

    confirmed: list[Pair] = []
    indeterminate: list[Pair] = []
    indet_alleles: set[str] = set()

    pair_count = 0
    for x, y in combinations_with_replacement(sorted(sigs), 2):
        pair_count += 1
        cmp = _compare_pair(observed, sigs[x], sigs[y])
        if cmp is ColumnComparison.INCOMPATIBLE:
            continue
        pair = (x, y) if x <= y else (y, x)
        if cmp is ColumnComparison.COMPATIBLE:
            confirmed.append(pair)
        else:  # INDETERMINATE
            indeterminate.append(pair)
            indet_alleles.update((x, y))
    limits.check_candidate_pairs(pair_count)

    # reference_policy handling of indeterminate candidates (main §4.0)
    if indeterminate and reference_policy is ReferencePolicy.STRICT:
        raise IndeterminateReferenceError(
            f"strict reference_policy: incomplete reference leaves "
            f"{len(indeterminate)} candidate pair(s) undecidable at {model.locus}; "
            f"cannot return an exact set",
            indeterminate_alleles=tuple(sorted(indet_alleles)),
        )
    # permissive / sequence_supplied keep the indeterminate candidates (flagged)

    return LocusCandidateSet(
        locus=model.locus,
        confirmed=tuple(sorted(set(confirmed))),
        indeterminate=tuple(sorted(set(indeterminate))),
        database_version=model.database_version,
        indeterminate_alleles=tuple(sorted(indet_alleles)),
    )


def to_genotype_list(
    candidates: LocusCandidateSet,
    *,
    representation_policy: RepresentationPolicy = RepresentationPolicy.EXACT,
    limits: ResourceLimits = DEFAULT_LIMITS,
) -> GenotypeList:
    """Reduce a candidate set to a GL ``GenotypeList`` (main §3.3).

    Exact by default: a lossless ``/`` factorization when the set is an exact
    Cartesian product, else exact ``|`` enumeration. The set used is
    confirmed + (for permissive) indeterminate.
    """
    pairs = list(candidates.all_pairs)
    if not pairs:
        raise RepresentationError(
            f"no candidate genotype at {candidates.locus} (empty set)"
        )
    limits.check_alternatives(len(pairs))
    return pairs_to_genotype_list(pairs)
