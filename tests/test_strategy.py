"""Strategy schema & loader tests (spec §3, QA §10.10)."""

from __future__ import annotations

from pathlib import Path

import pytest

from tissuetypeforge.strategy import (
    AssayClass,
    Resolution,
    StrategyValidationError,
    load_strategy,
    load_strategy_file,
)

STRAT_DIR = Path(__file__).resolve().parents[1] / "strategies"


def test_load_named_strategies_from_disk() -> None:
    for name in [
        "first_field",
        "second_field",
        "second_field_bc_phased",
        "mixed_registry",
    ]:
        s = load_strategy_file(STRAT_DIR / f"{name}.toml")
        assert s.name == name


def test_reference_database_recorded() -> None:
    # QA §10.10: reference database versions are recorded.
    s = load_strategy_file(STRAT_DIR / "second_field.toml")
    assert s.reference_database == "IPD-IMGT/HLA 3.55.0"
    cfg = s.config_for("HLA-A")
    assert cfg is not None
    assert cfg.reference_database == "IPD-IMGT/HLA 3.55.0"


def test_per_locus_override() -> None:
    s = load_strategy_file(STRAT_DIR / "mixed_registry.toml")
    assert s.config_for("HLA-C").resolution is Resolution.FIRST_FIELD
    assert s.config_for("HLA-A").resolution is Resolution.SECOND_FIELD
    assert s.config_for("HLA-DQB1").typed is False


def test_phase_group_parsed() -> None:
    s = load_strategy_file(STRAT_DIR / "second_field_bc_phased.toml")
    pg = s.phase_group_for("HLA-B")
    assert pg is not None
    assert pg.locus_set == frozenset({"HLA-B", "HLA-C"})
    assert s.phase_group_for("HLA-A") is None


def test_short_and_full_locus_keys_match() -> None:
    s = load_strategy(
        {"name": "t", "loci": {"A": {"resolution": "first_field"}}}
    )
    assert s.config_for("HLA-A").resolution is Resolution.FIRST_FIELD
    assert s.config_for("A").resolution is Resolution.FIRST_FIELD


def test_default_assay_class_is_idealized() -> None:
    s = load_strategy({"name": "t", "loci": {"A": {}}})
    assert s.config_for("A").assay_class is AssayClass.IDEALIZED


def test_invalid_resolution_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy({"name": "t", "loci": {"A": {"resolution": "bogus"}}})


def test_unknown_locus_key_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy({"name": "t", "loci": {"A": {"nonsense": 1}}})


def test_probability_out_of_range_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy({"name": "t", "loci": {"A": {"failure_probability": 1.5}}})


def test_single_locus_phase_group_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy({"name": "t", "phase_groups": [{"loci": ["B"]}]})


def test_overlapping_phase_groups_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy(
            {
                "name": "t",
                "phase_groups": [
                    {"loci": ["A", "B"]},
                    {"loci": ["B", "C"]},
                ],
            }
        )


def test_missing_name_rejected() -> None:
    with pytest.raises(StrategyValidationError):
        load_strategy({"loci": {"A": {}}})


def test_missing_file_raises() -> None:
    with pytest.raises(FileNotFoundError):
        load_strategy_file(STRAT_DIR / "does_not_exist.toml")
