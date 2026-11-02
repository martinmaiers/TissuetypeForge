"""Reference-data registry (Decision 4: reference-data selection).

A strategy selects a versioned reference artifact *by name* (``reference_name``);
the registry resolves that name to a concrete :class:`ColumnModel` (sequence) or
:class:`SerologyModel` (serology). This keeps strategies declarative — they name
an artifact, they do not embed one — and keeps TissueTypeForge independent of any
particular IPD distribution: the caller registers whatever pinned artifacts the
experiment uses (built by a separate data-prep tool, or hand-authored in tests).

The registry is explicit and in-process; there is no implicit network access.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from tissuetypeforge.reference.model import ColumnModel
from tissuetypeforge.serology.model import SerologyModel


class ReferenceLookupError(KeyError):
    """Raised when a strategy names a reference artifact not in the registry."""


@dataclass
class ReferenceRegistry:
    """Maps ``reference_name`` -> reference artifact, keyed by kind + locus.

    Sequence column models and serology models share a name-spaced registry so a
    single strategy can select both a sequence artifact and a serology artifact.
    Lookups are by ``(reference_name, locus)`` because one named database holds a
    per-locus artifact.
    """

    column_models: dict[tuple[str, str], ColumnModel] = field(default_factory=dict)
    serology_models: dict[tuple[str, str], SerologyModel] = field(default_factory=dict)

    def register_column_model(self, name: str, model: ColumnModel) -> None:
        self.column_models[(name, model.locus)] = model

    def register_serology_model(self, name: str, model: SerologyModel) -> None:
        self.serology_models[(name, model.locus)] = model

    def column_model(self, name: str, locus: str) -> ColumnModel:
        try:
            return self.column_models[(name, locus)]
        except KeyError as exc:
            raise ReferenceLookupError(
                f"no sequence reference {name!r} for locus {locus!r}; "
                f"register one via ReferenceRegistry.register_column_model"
            ) from exc

    def serology_model(self, name: str, locus: str) -> SerologyModel:
        try:
            return self.serology_models[(name, locus)]
        except KeyError as exc:
            raise ReferenceLookupError(
                f"no serology reference {name!r} for locus {locus!r}; "
                f"register one via ReferenceRegistry.register_serology_model"
            ) from exc
