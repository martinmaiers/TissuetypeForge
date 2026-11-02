"""M2B — policies, resource limits, and error types for the sequence engine.

Implements the approved final decisions:

* reference policy strict(default)/permissive/sequence_supplied, where ``strict``
  is fail-closed: it raises rather than silently excluding indeterminate
  candidates and calling the result exact (main §4.0).
* representation policy exact(default)/superset/annotate_only (main §3.3).
* configurable resource limits producing explicit ``ResourceLimitError``
  (main §3.3.3, D4).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class ReferencePolicy(StrEnum):
    STRICT = "strict"
    PERMISSIVE = "permissive"
    SEQUENCE_SUPPLIED = "sequence_supplied"


class RepresentationPolicy(StrEnum):
    EXACT = "exact"
    SUPERSET = "superset"
    ANNOTATE_ONLY = "annotate_only"


class ObservationError(ValueError):
    """Base class for sequence-observation errors."""


class IndeterminateReferenceError(ObservationError):
    """Raised by ``strict`` policy when incomplete reference prevents an exact set.

    Carries the alleles whose membership could not be decided, so the caller can
    see exactly why exactness failed (never a silently reduced set).
    """

    def __init__(self, message: str, *, indeterminate_alleles: tuple[str, ...]) -> None:
        super().__init__(message)
        self.indeterminate_alleles = indeterminate_alleles


class MissingReferenceError(ObservationError):
    """Raised when the TRUE allele lacks reference over the interrogated region."""


class RepresentationError(ObservationError):
    """Raised under ``exact`` policy when a result cannot be expressed in GL."""


class ResourceLimitError(ObservationError):
    """Raised when a configured resource limit is exceeded (D4)."""

    def __init__(self, message: str, *, kind: str) -> None:
        super().__init__(message)
        self.kind = kind


@dataclass(frozen=True)
class ResourceLimits:
    """Configurable bounds; exceeding one raises :class:`ResourceLimitError`.

    ``None`` means unbounded for that dimension. Defaults are generous but
    finite so a pathological input fails explicitly rather than hanging.
    """

    max_universe: int | None = 100_000
    max_candidate_pairs: int | None = 1_000_000
    max_alternatives: int | None = 10_000
    max_joint_molecules: int | None = 1_000_000

    def check_universe(self, n: int) -> None:
        if self.max_universe is not None and n > self.max_universe:
            raise ResourceLimitError(
                f"candidate universe size {n} exceeds max_universe {self.max_universe}",
                kind="universe",
            )

    def check_candidate_pairs(self, n: int) -> None:
        if self.max_candidate_pairs is not None and n > self.max_candidate_pairs:
            raise ResourceLimitError(
                f"candidate pairs {n} exceeds max_candidate_pairs "
                f"{self.max_candidate_pairs}",
                kind="candidate_pairs",
            )

    def check_alternatives(self, n: int) -> None:
        if self.max_alternatives is not None and n > self.max_alternatives:
            raise ResourceLimitError(
                f"genotype alternatives {n} exceeds max_alternatives "
                f"{self.max_alternatives}",
                kind="alternatives",
            )

    def check_joint_molecules(self, n: int) -> None:
        if self.max_joint_molecules is not None and n > self.max_joint_molecules:
            raise ResourceLimitError(
                f"joint molecule pairings {n} exceeds max_joint_molecules "
                f"{self.max_joint_molecules}",
                kind="joint_molecules",
            )


#: Shared default limits singleton (avoids a function call in argument defaults).
DEFAULT_LIMITS = ResourceLimits()
