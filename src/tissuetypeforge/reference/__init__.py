"""Versioned reference artifacts and sequence signatures (M2A)."""

from tissuetypeforge.reference.model import Base, ColumnModel, FeatureKey
from tissuetypeforge.reference.registry import ReferenceLookupError, ReferenceRegistry
from tissuetypeforge.reference.signature import (
    ColumnComparison,
    compare_column,
    interrogated_columns,
    signature,
)

__all__ = [
    "Base",
    "ColumnComparison",
    "ColumnModel",
    "FeatureKey",
    "ReferenceLookupError",
    "ReferenceRegistry",
    "compare_column",
    "interrogated_columns",
    "signature",
]
