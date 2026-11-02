"""Serology observation backend (M2C)."""

from tissuetypeforge.serology.backend import SerologyObservation, observe_serology
from tissuetypeforge.serology.model import (
    Assignment,
    BroadSplitLattice,
    Confidence,
    Epitope,
    ResolutionLevel,
    SerologicAntigenReport,
    SerologyModel,
    Specificity,
)

__all__ = [
    "Assignment",
    "BroadSplitLattice",
    "Confidence",
    "Epitope",
    "ResolutionLevel",
    "SerologicAntigenReport",
    "SerologyModel",
    "SerologyObservation",
    "Specificity",
    "observe_serology",
]
