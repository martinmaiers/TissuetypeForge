"""Load typing strategies from TOML (spec §3).

TOML is chosen over YAML to avoid an external dependency: ``tomllib`` is in the
standard library from Python 3.11. The schema:

    name = "second_field"
    description = "..."
    reference_database = "IPD-IMGT/HLA 3.55.0"

    [defaults]
    typed = true
    assay_class = "idealized"
    resolution = "second_field"
    failure_probability = 0.0
    error_probability = 0.0

    [loci.A]
    resolution = "first_field"

    [[phase_groups]]
    loci = ["B", "C"]
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

from tissuetypeforge.strategy.schema import (
    AssayClass,
    LocusConfig,
    PhaseGroup,
    Resolution,
    Strategy,
    StrategyValidationError,
)

_KNOWN_LOCUS_KEYS = {
    "typed",
    "assay_class",
    "resolution",
    "reference_database",
    "reference_name",
    "interrogated_features",
    "reference_policy",
    "serology_level",
    "observed_epitopes",
    "sequence_regions",
    "ambiguity_rules",
    "failure_probability",
    "error_probability",
    "assay_params",
}


def _coerce_enum(enum_cls: type, value: Any, field_name: str) -> Any:
    try:
        return enum_cls(value)
    except ValueError as exc:
        allowed = ", ".join(e.value for e in enum_cls)
        raise StrategyValidationError(
            f"Invalid {field_name} {value!r}; allowed: {allowed}"
        ) from exc


def _build_locus_config(
    locus: str, raw: dict[str, Any], defaults: dict[str, Any]
) -> LocusConfig:
    merged: dict[str, Any] = {**defaults, **raw}
    unknown = set(raw) - _KNOWN_LOCUS_KEYS
    if unknown:
        raise StrategyValidationError(
            f"Unknown keys {sorted(unknown)} in config for locus {locus}"
        )
    assay_class = _coerce_enum(
        AssayClass, merged.get("assay_class", "idealized"), "assay_class"
    )
    resolution = _coerce_enum(
        Resolution, merged.get("resolution", "second_field"), "resolution"
    )
    return LocusConfig(
        locus=locus,
        typed=bool(merged.get("typed", True)),
        assay_class=assay_class,
        resolution=resolution,
        reference_database=merged.get("reference_database"),
        reference_name=merged.get("reference_name"),
        interrogated_features=tuple(merged.get("interrogated_features", ()) or ()),
        reference_policy=merged.get("reference_policy"),
        serology_level=merged.get("serology_level", "split"),
        observed_epitopes=tuple(merged.get("observed_epitopes", ()) or ()),
        sequence_regions=tuple(merged.get("sequence_regions", ()) or ()),
        ambiguity_rules=merged.get("ambiguity_rules"),
        failure_probability=float(merged.get("failure_probability", 0.0)),
        error_probability=float(merged.get("error_probability", 0.0)),
        assay_params=dict(merged.get("assay_params", {}) or {}),
    )


def load_strategy(data: dict[str, Any]) -> Strategy:
    """Build a :class:`Strategy` from an already-parsed TOML mapping."""
    name = data.get("name")
    if not name:
        raise StrategyValidationError("Strategy must have a 'name'")

    defaults_raw: dict[str, Any] = dict(data.get("defaults", {}))
    unknown_defaults = set(defaults_raw) - _KNOWN_LOCUS_KEYS
    if unknown_defaults:
        raise StrategyValidationError(
            f"Unknown keys {sorted(unknown_defaults)} in [defaults]"
        )

    reference_database = data.get("reference_database")
    # Thread a top-level reference_database into defaults if not overridden.
    if reference_database is not None:
        defaults_raw.setdefault("reference_database", reference_database)

    loci_raw: dict[str, Any] = dict(data.get("loci", {}))
    loci: dict[str, LocusConfig] = {}
    for locus_key, raw in loci_raw.items():
        if not isinstance(raw, dict):
            raise StrategyValidationError(
                f"[loci.{locus_key}] must be a table, got {type(raw).__name__}"
            )
        loci[locus_key] = _build_locus_config(locus_key, raw, defaults_raw)

    # Default config applied to loci present in truth but absent from [loci].
    default_locus = _build_locus_config("<default>", {}, defaults_raw)

    phase_groups_raw = data.get("phase_groups", [])
    if not isinstance(phase_groups_raw, list):
        raise StrategyValidationError("'phase_groups' must be an array of tables")
    phase_groups: list[PhaseGroup] = []
    for pg in phase_groups_raw:
        if not isinstance(pg, dict) or "loci" not in pg:
            raise StrategyValidationError(
                "Each phase group must be a table with a 'loci' array"
            )
        phase_groups.append(PhaseGroup(tuple(pg["loci"])))

    return Strategy(
        name=name,
        loci=loci,
        phase_groups=tuple(phase_groups),
        reference_database=reference_database,
        description=data.get("description", ""),
        default_locus=default_locus,
        reference_policy=data.get("reference_policy", "strict"),
        resource_limits=dict(data.get("resource_limits", {}) or {}),
    )


def load_strategy_file(path: str | Path) -> Strategy:
    """Load and validate a strategy from a TOML file."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Strategy file not found: {p}")
    with p.open("rb") as fh:
        data = tomllib.load(fh)
    return load_strategy(data)
