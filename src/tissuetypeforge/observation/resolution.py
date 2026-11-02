"""Idealized resolution observations (spec §5, M1 scope).

These model what the *reported name* would be under an idealized assumption:
that the method resolves exactly to a nomenclature field with no measurement
ambiguity. This is a baseline convenience model, NOT an empirical assay (spec
§5, §12). Sequence-derived ambiguity, serology, intermediate resolution, and
full-gene resolution are deferred to M2 and raise ``NotImplementedError`` here
rather than being faked.

Idealized modes do not introduce ambiguity: a field-truncated allele is
reported as a single allele, not a MAC/ambiguity set, because the idealized
model asserts exact resolution to that field (README scientific assumptions).
"""

from __future__ import annotations

from tissuetypeforge.glstring import Allele, AlleleList
from tissuetypeforge.nomenclature import truncate_to_fields
from tissuetypeforge.strategy.schema import Resolution

_M2_MESSAGE = (
    "Resolution mode {mode!r} requires assay-derived modeling (reference "
    "sequence or serology mappings) and is deferred to Milestone 2. M1 "
    "implements only idealized 'allele', 'first_field', and 'second_field'. "
    "Refusing to fake an empirical transformation (spec §12)."
)


def resolve_allele_list(seg: AlleleList, resolution: Resolution) -> AlleleList:
    """Apply an idealized resolution to every allele in an ambiguity segment.

    Truncation can collapse distinct alleles to the same field-level name; the
    result is de-duplicated while preserving order so an idealized observation
    does not report a spurious ``a/a`` ambiguity.
    """
    if not resolution.is_idealized_m1:
        raise NotImplementedError(_M2_MESSAGE.format(mode=resolution.value))

    resolved: list[Allele] = [_resolve_one(a, resolution) for a in seg.alleles]
    deduped: list[Allele] = []
    seen: set[str] = set()
    for a in resolved:
        if a.name not in seen:
            seen.add(a.name)
            deduped.append(a)
    return AlleleList(tuple(deduped))


def _resolve_one(allele: Allele, resolution: Resolution) -> Allele:
    if resolution is Resolution.ALLELE:
        return allele
    if resolution is Resolution.FIRST_FIELD:
        return truncate_to_fields(allele, 1)
    if resolution is Resolution.SECOND_FIELD:
        return truncate_to_fields(allele, 2)
    raise NotImplementedError(_M2_MESSAGE.format(mode=resolution.value))
