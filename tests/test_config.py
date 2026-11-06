"""Tests for the TOML/dict population config loader (REQ 6, deliverables)."""

from __future__ import annotations

from pathlib import Path

import pytest

from tissuetypeforge import load_population, load_population_file
from tissuetypeforge.population import ProfileValidationError

POP_DIR = Path(__file__).resolve().parents[1] / "populations"


def test_load_example_population_file() -> None:
    pop = load_population_file(POP_DIR / "registry_example.toml")
    assert pop.name == "registry-example"
    assert pop.loci == ("A", "C", "B", "DRB1")
    assert [p.id for p in pop.profiles] == ["legacy", "intermediate", "modern"]
    # 0-field for C in the legacy profile is preserved.
    assert pop.profiles[0].mode_for("C").value == "0-field"


def test_missing_population_file_raises() -> None:
    with pytest.raises(FileNotFoundError):
        load_population_file(POP_DIR / "does_not_exist.toml")


def test_load_population_rejects_non_mapping() -> None:
    with pytest.raises(ProfileValidationError):
        load_population([])  # type: ignore[arg-type]


def test_load_population_from_dict_roundtrip() -> None:
    data = {
        "name": "x",
        "loci": ["A", "B"],
        "profiles": [
            {"id": "only", "weight": 1.0, "resolution": {"A": "2-field", "B": "2-field"}},
        ],
    }
    pop = load_population(data)
    assert pop.name == "x"
    assert pop.normalized_weights() == [1.0]
