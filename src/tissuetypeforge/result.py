"""Structured output of an observation (spec §2, §7).

The output contract is designed so the Milestone 3 error layer slots in without
changing shape. It carries, side by side:

    true_diplotype          the known truth GL String (retained so mistyping is
                            measurable — spec §7)
    expected_observation    the error-free observation under the assay
    observed                the actual simulated observation (== expected in M1)
    error_category          none / information_loss / missing / incorrect
    error_events            per-event provenance (empty in M1)
    metadata                how the observation was produced
    loss_record             what information the assay could not observe

In Milestone 1 there is no stochastic error, so ``observed`` equals
``expected_observation`` and ``error_category`` is always ``NONE``. Ordinary
assay ambiguity and information loss are recorded in ``loss_record`` and are
*not* labeled as mistyping (spec §7).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class ErrorCategory(StrEnum):
    """Mutually distinct outcome categories (spec §7).

    Information loss / ambiguity, missing typing, and incorrect typing are not
    equivalent and must be distinguished.
    """

    NONE = "none"
    INFORMATION_LOSS = "information_loss"
    MISSING = "missing"
    INCORRECT = "incorrect"


@dataclass(frozen=True)
class LocusObservation:
    """Per-locus record of what happened during observation."""

    locus: str
    typed: bool
    resolution: str | None
    assay_class: str | None
    reference_database: str | None
    #: Human-readable notes on information retained / lost at this locus.
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ErrorEvent:
    """Provenance for a single simulated error (populated from M3)."""

    locus: str
    mechanism: str
    detail: str


@dataclass(frozen=True)
class ObservationMetadata:
    """Describes how the observation was produced (spec §2)."""

    strategy_name: str
    reference_databases: tuple[str, ...]
    typed_loci: tuple[str, ...]
    missing_loci: tuple[str, ...]
    retained_phase_groups: tuple[tuple[str, ...], ...]
    locus_observations: tuple[LocusObservation, ...]
    seed: int | None = None


@dataclass(frozen=True)
class ObservationResult:
    """Full structured result of observing one truth genotype."""

    true_diplotype: str
    expected_observation: str
    observed: str
    error_category: ErrorCategory
    metadata: ObservationMetadata
    error_events: tuple[ErrorEvent, ...] = ()
    loss_record: tuple[str, ...] = field(default_factory=tuple)
    #: Serologic antigen reports, one per serology-typed locus (M2, Decision 4).
    #: Kept SEPARATE from the observed GL String: antigen tokens are not GL
    #: alleles. Empty for non-serologic observations (backward-compatible).
    antigen_reports: tuple[object, ...] = ()

    @property
    def had_error(self) -> bool:
        return self.error_category is not ErrorCategory.NONE
