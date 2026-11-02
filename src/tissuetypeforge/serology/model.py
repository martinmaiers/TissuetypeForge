"""M2C — versioned serology reference model (serology companion design §2-§5).

Serology is NOT a truncated DNA name. This module holds the versioned
allele->specificity assignment, confidence, and the broad/split lattice. Antigen
and epitope tokens are a *distinct type* from nomenclature alleles and must never
be serialized into a GL String (companion §7).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class Confidence(StrEnum):
    """Assignment confidence (companion §2, from rel_dna_ser semantics)."""

    ASSIGNED = "assigned"          # '+'  expressed, confirmed
    POSSIBLE = "possible"          # '?'  expected but unconfirmed
    UNKNOWN = "unknown"            # blank cell: no assignment on record
    NOT_EXPRESSED = "not_expressed"  # '0'  null / blank antigen


@dataclass(frozen=True)
class Specificity:
    """A serologic antigen token, e.g. ``B7``, ``B62``. NOT a GL allele."""

    name: str


@dataclass(frozen=True)
class Epitope:
    """A public-epitope token, e.g. ``Bw4``, ``Bw6``. Separate axis; NOT a GL allele."""

    name: str


class ResolutionLevel(StrEnum):
    BROAD = "broad"
    SPLIT = "split"


@dataclass(frozen=True)
class Assignment:
    """Serologic assignment for one allele."""

    specificities: frozenset[Specificity]
    confidence: Confidence
    epitopes: frozenset[Epitope] = frozenset()


@dataclass(frozen=True)
class BroadSplitLattice:
    """Broad/split relationships (companion §3)."""

    broad_of: dict[Specificity, Specificity]  # split -> broad (or self)

    def project(self, spec: Specificity, level: ResolutionLevel) -> Specificity:
        if level is ResolutionLevel.SPLIT:
            return spec
        return self.broad_of.get(spec, spec)


@dataclass(frozen=True)
class SerologyModel:
    """Versioned serology artifact (companion §2)."""

    locus: str
    database_version: str
    assignment: dict[str, Assignment]
    lattice: BroadSplitLattice = field(
        default_factory=lambda: BroadSplitLattice(broad_of={})
    )

    def alleles(self) -> tuple[str, ...]:
        return tuple(sorted(self.assignment))

    def get(self, allele: str) -> Assignment | None:
        return self.assignment.get(allele)

    def projected_specificities(
        self, allele: str, level: ResolutionLevel
    ) -> frozenset[Specificity] | None:
        """Antigen set an allele presents at ``level``; ``None`` if allele absent.

        A ``NOT_EXPRESSED`` (null) allele yields the empty set (blank). An
        ``UNKNOWN`` allele returns ``None`` is NOT used here — absence from the
        model is ``None``; UNKNOWN confidence is a *present* assignment with no
        asserted antigen, handled by the backend as indeterminate.
        """
        a = self.assignment.get(allele)
        if a is None:
            return None
        if a.confidence is Confidence.NOT_EXPRESSED:
            return frozenset()
        return frozenset(self.lattice.project(s, level) for s in a.specificities)


@dataclass(frozen=True)
class SerologicAntigenReport:
    """The observed antigen report — a structured object, NOT a GL String.

    Antigen/epitope tokens live here only; they are never placed in the
    compatible-DNA GL String (companion §7).
    """

    locus: str
    specificities: frozenset[Specificity]
    resolution_level: ResolutionLevel
    public_epitopes: frozenset[Epitope]
    database_version: str
    uses_possible_assignment: bool = False
    uses_unknown_assignment: bool = False
    blank_present: bool = False
