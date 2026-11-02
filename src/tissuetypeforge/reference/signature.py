"""M2A — observable sequence signatures and three-valued comparison.

A *signature* is the projection of an allele's reference row onto the columns an
assay interrogates (main design §2.2). Comparison of an observed base-set
against a candidate base is three-valued (main §4.1, D2): a reference gap is
INDETERMINATE and can never confirm compatibility.
"""

from __future__ import annotations

from enum import StrEnum

from tissuetypeforge.reference.model import Base, ColumnModel, FeatureKey


class ColumnComparison(StrEnum):
    """Three-valued outcome of comparing an observed column to a candidate."""

    COMPATIBLE = "compatible"
    INCOMPATIBLE = "incompatible"
    INDETERMINATE = "indeterminate"


def signature(
    model: ColumnModel, allele: str, columns: tuple[int, ...]
) -> tuple[Base, ...] | None:
    """Project ``allele``'s reference row onto ``columns``.

    Returns ``None`` if the allele is absent from the reference database (a
    distinct condition from a signature containing INDETERMINATE — that means
    present-but-unsequenced; main §4). The caller decides how absence is handled
    per ``reference_policy``.
    """
    row = model.row(allele)
    if row is None:
        return None
    return tuple(row[i] for i in columns)


def compare_column(observed: frozenset[Base], candidate_base: Base) -> ColumnComparison:
    """Compare an observed base-set at one column to a candidate's base (§4.1).

    ``observed`` is the set of bases the diploid assay saw at this column across
    both chromosomes (cardinality 1 = homozygous, 2 = heterozygous; §2.3.1).
    A candidate INDETERMINATE base (reference gap / missing data) is never a
    confirmed match.
    """
    if candidate_base is Base.INDETERMINATE:
        return ColumnComparison.INDETERMINATE
    if Base.INDETERMINATE in observed:
        # Truth side itself unknown at this interrogated column.
        return ColumnComparison.INDETERMINATE
    if candidate_base in observed:
        return ColumnComparison.COMPATIBLE
    return ColumnComparison.INCOMPATIBLE


def interrogated_columns(
    model: ColumnModel, features: frozenset[FeatureKey]
) -> tuple[int, ...]:
    """Columns the assay reads = union of the interrogated features' columns."""
    cols = model.columns_for_features(features)
    if not cols:
        raise ValueError(
            f"No columns found for interrogated features {sorted(map(str, features))} "
            f"in locus {model.locus}"
        )
    return cols
