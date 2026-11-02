"""Shared GL-reduction helpers (exact ``|`` / lossless ``/``) — main §3.3.

Extracted into a dependency-free module so both the sequence and serology
backends can reuse it without importing each other (avoids a circular import).
"""

from __future__ import annotations

from itertools import product

from tissuetypeforge.glstring import Allele, AlleleList, Genotype, GenotypeList, PhasedGenotype

Pair = tuple[str, str]


def factor_cartesian(pairs: tuple[Pair, ...]) -> tuple[list[str], list[str]] | None:
    """If ``pairs`` equals an unordered Cartesian product LxR, return (L, R).

    Only then is a ``/``-factorization lossless (main §3.3 step 2); otherwise
    return ``None`` and the caller emits exact ``|``. Bounded ``3^n`` bipartition
    search over the involved alleles (Decision 2: bounded; fall back to ``|``).
    """
    if not pairs:
        return None
    alleles = sorted({a for p in pairs for a in p})
    if len(alleles) > 12:
        return None
    target = set(pairs)

    def product_pairs(left: list[str], right: list[str]) -> set[Pair]:
        return {(a, b) if a <= b else (b, a) for a in left for b in right}

    for assign in product("LRB", repeat=len(alleles)):
        left = [alleles[i] for i, s in enumerate(assign) if s in ("L", "B")]
        right = [alleles[i] for i, s in enumerate(assign) if s in ("R", "B")]
        if not left or not right:
            continue
        if product_pairs(left, right) == target:
            return sorted(set(left)), sorted(set(right))
    return None


def pair_to_genotype(pair: Pair) -> Genotype:
    x, y = pair
    return Genotype(
        (
            PhasedGenotype((AlleleList((Allele(x),)),)),
            PhasedGenotype((AlleleList((Allele(y),)),)),
        )
    )


def pairs_to_genotype_list(pairs: list[Pair]) -> GenotypeList:
    """Reduce an unordered-pair set to a GenotypeList: lossless ``/`` else ``|``."""
    if not pairs:
        raise ValueError("empty candidate pair set")
    if len(pairs) == 1:
        return GenotypeList((pair_to_genotype(pairs[0]),))
    factored = factor_cartesian(tuple(pairs))
    if factored is not None:
        left, right = factored
        geno = Genotype(
            (
                PhasedGenotype((AlleleList(tuple(Allele(a) for a in left)),)),
                PhasedGenotype((AlleleList(tuple(Allele(a) for a in right)),)),
            )
        )
        return GenotypeList((geno,))
    return GenotypeList(tuple(pair_to_genotype(p) for p in pairs))
