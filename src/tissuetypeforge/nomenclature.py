"""Nomenclature helpers for the SIP ``glstring`` ``Allele`` type.

The library's :class:`glstring.Allele` is intentionally minimal (``name``,
``locus``, ``fields``) and does not provide field truncation. TissueTypeForge's
idealized resolution modes need truncation, so it lives here as a free function
with the scientifically correct expression-suffix handling carried over from
Milestone 1.
"""

from __future__ import annotations

from tissuetypeforge.glstring import Allele


def truncate_to_fields(allele: Allele, n: int) -> Allele:
    """Return a new :class:`Allele` truncated to the first ``n`` fields.

    The WHO expression suffix (N, L, S, C, A, Q) is written on the last field of
    the full name. Truncation keeps the first ``n`` fields exactly as written: a
    surviving suffix-bearing field keeps its suffix; a truncated-away one drops
    it. The suffix is never migrated onto a different field (asserting e.g. null
    expression at a resolution that did not observe it would be scientifically
    wrong). If there are fewer than ``n`` fields, the allele is returned
    unchanged.
    """
    if n < 1:
        raise ValueError("Field count must be >= 1")
    flds = allele.fields
    if not flds or len(flds) <= n:
        return allele
    kept = flds[:n]
    return Allele(f"{allele.locus}*{':'.join(kept)}")
