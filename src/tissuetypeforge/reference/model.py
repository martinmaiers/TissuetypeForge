"""M2A — versioned reference artifacts (column model).

Implements decision D1 (pinned, checksummed per-locus column model) and D2
(``gap`` and ``indeterminate`` are distinct symbols). This is the reference
substrate the sequence backend reads; it is deliberately independent of any live
IPD-IMGT/HLA file format so tests can use small hand-authored artifacts.

Three information layers (main design §4) are kept distinct here:

* biological truth — not represented in the artifact (it is the input genotype);
* **available reference knowledge** — THIS module: per-allele bases over a common
  alignment frame, each column flagged PRESENT / INDETERMINATE;
* observable information — derived by projecting onto interrogated columns
  (:mod:`tissuetypeforge.reference.signature`).

A column holds one :class:`Base` symbol. Critically (D2):

* ``GAP`` is a *sequenced* alignment gap (a real, observed deletion) — evidence
  that can distinguish alleles.
* ``INDETERMINATE`` is "no reference base here" (missing data) — never evidence
  of identity.

The column-model builder (out of scope: a data-prep tool) is responsible for the
invariant that a gap arising from missing data is encoded ``INDETERMINATE``, not
``GAP``.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from enum import StrEnum


class Base(StrEnum):
    """A single alignment-column symbol (D2)."""

    A = "A"
    C = "C"
    G = "G"
    T = "T"
    #: A sequenced deletion (observed absence of base) — discriminating evidence.
    GAP = "-"
    #: No reference base on record for this allele/column — missing data.
    INDETERMINATE = "?"


@dataclass(frozen=True)
class FeatureKey:
    """Identifies a sequence feature, e.g. ``("exon", 2)`` or ``("intron", 1)``."""

    kind: str
    index: int

    def __str__(self) -> str:
        return f"{self.kind}{self.index}"


@dataclass(frozen=True)
class ColumnModel:
    """A pinned, versioned per-locus alignment-column model (D1).

    ``rows`` maps an allele name to a tuple of :class:`Base`, one per alignment
    column, all of equal length ``n_columns``. ``feature_of_column`` maps each
    column index to its :class:`FeatureKey`, so an assay can select the columns
    of the features it interrogates.
    """

    locus: str
    database_version: str
    feature_of_column: tuple[FeatureKey, ...]
    rows: dict[str, tuple[Base, ...]]
    checksum: str = field(default="")

    def __post_init__(self) -> None:
        n = len(self.feature_of_column)
        for allele, row in self.rows.items():
            if len(row) != n:
                raise ValueError(
                    f"Allele {allele} row length {len(row)} != n_columns {n}"
                )
        if not self.checksum:
            object.__setattr__(self, "checksum", self._compute_checksum())

    @property
    def n_columns(self) -> int:
        return len(self.feature_of_column)

    def columns_for_features(self, features: frozenset[FeatureKey]) -> tuple[int, ...]:
        """Return the ordered column indices belonging to any of ``features``."""
        return tuple(
            i for i, fk in enumerate(self.feature_of_column) if fk in features
        )

    def alleles(self) -> tuple[str, ...]:
        """All allele names in the model, in sorted order (deterministic)."""
        return tuple(sorted(self.rows))

    def row(self, allele: str) -> tuple[Base, ...] | None:
        """Return the base row for ``allele``, or ``None`` if not in the model.

        ``None`` means the allele is *absent from the reference database* — a
        distinct condition from a row full of INDETERMINATE (present in the
        database but unsequenced). Callers must not conflate the two.
        """
        return self.rows.get(allele)

    def _compute_checksum(self) -> str:
        h = hashlib.sha256()
        h.update(self.locus.encode())
        h.update(self.database_version.encode())
        h.update(repr(self.feature_of_column).encode())
        for allele in sorted(self.rows):
            h.update(allele.encode())
            h.update("".join(self.rows[allele]).encode())
        return h.hexdigest()
