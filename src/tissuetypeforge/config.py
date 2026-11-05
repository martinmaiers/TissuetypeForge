"""Load a :class:`PopulationTypingModel` from TOML or an in-memory mapping.

TOML is used (stdlib ``tomllib``, Python 3.11+) so no new dependency is needed.
The mapping shape mirrors the requirement's example exactly::

    name = "registry-x"
    loci = ["A", "C", "B", "DRB1"]   # optional; order is preserved

    [[profiles]]
    id = "legacy"
    weight = 0.20
    resolution = { A = "1-field", B = "1-field", C = "0-field", DRB1 = "1-field" }

    [[profiles]]
    id = "modern"
    weight = 0.80
    resolution = { A = "2-field", B = "2-field", C = "2-field", DRB1 = "2-field" }

The same structure as a Python ``dict`` (e.g. parsed from JSON) is accepted by
:func:`load_population`.
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

from tissuetypeforge.population import (
    PopulationTypingModel,
    ProfileValidationError,
    build_population,
)


def load_population(data: dict[str, Any]) -> PopulationTypingModel:
    """Build and validate a population model from a parsed mapping (REQ 5)."""
    if not isinstance(data, dict):
        raise ProfileValidationError("population config must be a mapping/table")
    return build_population(data)


def load_population_file(path: str | Path) -> PopulationTypingModel:
    """Load and validate a population model from a TOML file."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Population config file not found: {p}")
    with p.open("rb") as fh:
        data = tomllib.load(fh)
    return load_population(data)
