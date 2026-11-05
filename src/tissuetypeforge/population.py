"""Population-level joint typing-profile distributions (REQ 2, 3, 5).

A **typing profile** is a complete, named assignment of a resolution mode to
*every* configured locus. A **population typing model** is a categorical
distribution over profiles: each simulated donor draws exactly one profile, and
the resolution choices across loci are therefore *joint*, never sampled
independently per locus.

This models how real donor registries contain groups of donors typed under
different historical protocols — a donor typed by a "legacy" lab has a
*correlated* resolution pattern across all loci, not an independent coin-flip at
each locus.

Separation of concerns (REQ 3)
------------------------------
The population's *genetic composition* (which true genotypes occur, and how
often) is a completely separate concept from its *typing-profile distribution*
(how donors were typed). This module models only the latter. No association
between genotype and typing profile is introduced unless a caller explicitly
supplies one (none is supported here by design — the sampler is genotype-blind).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from tissuetypeforge.resolution import ResolutionError, ResolutionMode, parse_resolution

#: Tolerance for the "weights sum to 1" check (REQ 5).
_WEIGHT_SUM_TOLERANCE = 1e-6


class ProfileValidationError(ValueError):
    """Raised when a profile or population configuration is invalid (REQ 5)."""


def _normalize_locus(locus: str) -> str:
    """Canonicalize a locus key to its bare form (``HLA-A`` -> ``A``).

    Profiles may be written with either the bare (``A``) or prefixed (``HLA-A``)
    form; both refer to the same locus. Observation matches on the bare form.
    """
    return locus.split("-", 1)[-1]


@dataclass(frozen=True)
class TypingProfile:
    """A complete per-locus resolution assignment with an id and weight.

    * ``id`` — unique within a population.
    * ``weight`` — unnormalized categorical weight (finite, >= 0). The sampler
      normalizes weights across the population.
    * ``resolution`` — maps every configured locus (bare form) to a
      :class:`ResolutionMode`.
    """

    id: str
    weight: float
    resolution: dict[str, ResolutionMode]

    def mode_for(self, locus: str) -> ResolutionMode | None:
        """Resolution mode for ``locus`` (accepts bare or ``HLA-`` form)."""
        return self.resolution.get(_normalize_locus(locus))


@dataclass(frozen=True)
class PopulationTypingModel:
    """A categorical distribution of typing profiles for one population/cohort.

    ``loci`` is the ordered tuple of configured loci (bare form). Every profile
    must assign a resolution mode to exactly these loci (REQ 5). ``name``
    identifies the cohort so different registries can carry different
    distributions (REQ 3).
    """

    name: str
    loci: tuple[str, ...]
    profiles: tuple[TypingProfile, ...]

    @property
    def total_weight(self) -> float:
        return math.fsum(p.weight for p in self.profiles)

    def normalized_weights(self) -> list[float]:
        """Return per-profile probabilities normalized to sum to 1."""
        total = self.total_weight
        return [p.weight / total for p in self.profiles]


def _coerce_mode(value: object, *, profile_id: str, locus: str) -> ResolutionMode:
    try:
        return parse_resolution(value)
    except ResolutionError as exc:
        raise ProfileValidationError(
            f"profile {profile_id!r}, locus {locus!r}: {exc}"
        ) from exc


def build_profile(raw: dict, *, loci: tuple[str, ...]) -> TypingProfile:
    """Build and validate a single :class:`TypingProfile` against ``loci``.

    Validates (REQ 5): the id is present and non-empty; the weight is finite and
    nonnegative; the resolution covers exactly the configured loci with values
    from the controlled vocabulary; unknown loci are rejected with a clear error.
    """
    if "id" not in raw or not str(raw["id"]).strip():
        raise ProfileValidationError("each profile must have a non-empty 'id'")
    pid = str(raw["id"])

    if "weight" not in raw:
        raise ProfileValidationError(f"profile {pid!r} must specify a 'weight'")
    try:
        weight = float(raw["weight"])
    except (TypeError, ValueError) as exc:
        raise ProfileValidationError(
            f"profile {pid!r}: weight {raw['weight']!r} is not a number"
        ) from exc
    if not math.isfinite(weight):
        raise ProfileValidationError(
            f"profile {pid!r}: weight must be finite, got {weight}"
        )
    if weight < 0:
        raise ProfileValidationError(
            f"profile {pid!r}: weight must be nonnegative, got {weight}"
        )

    res_raw = raw.get("resolution")
    if not isinstance(res_raw, dict):
        raise ProfileValidationError(
            f"profile {pid!r}: 'resolution' must be a table mapping locus -> mode"
        )

    configured = set(loci)
    provided = {_normalize_locus(k) for k in res_raw}

    unknown = provided - configured
    if unknown:
        raise ProfileValidationError(
            f"profile {pid!r}: resolution names unknown loci {sorted(unknown)}; "
            f"configured loci are {sorted(configured)}"
        )
    missing = configured - provided
    if missing:
        raise ProfileValidationError(
            f"profile {pid!r}: resolution is missing configured loci "
            f"{sorted(missing)}; every profile must specify every configured locus"
        )

    resolution = {
        _normalize_locus(k): _coerce_mode(v, profile_id=pid, locus=k)
        for k, v in res_raw.items()
    }
    return TypingProfile(id=pid, weight=weight, resolution=resolution)


def build_population(data: dict) -> PopulationTypingModel:
    """Build and validate a :class:`PopulationTypingModel` from a mapping (REQ 5).

    Expected shape (matching the requirement example)::

        {
          "name": "cohort-x",             # optional, defaults to "population"
          "loci": ["A", "C", "B", "DRB1"],# optional; inferred from profiles if
                                          # omitted (must then agree across all)
          "profiles": [
            {"id": "legacy", "weight": 0.2,
             "resolution": {"A": "1-field", "B": "1-field",
                            "C": "0-field", "DRB1": "1-field"}},
            ...
          ]
        }

    Validation (REQ 5):

    * profile ids are unique;
    * weights finite, nonnegative;
    * at least one profile has positive weight;
    * weights sum to 1 within ``_WEIGHT_SUM_TOLERANCE``;
    * every profile specifies every configured locus;
    * resolution values belong to the controlled vocabulary;
    * unknown loci / unsupported modes fail with a clear message.
    """
    name = str(data.get("name", "population"))
    profiles_raw = data.get("profiles")
    if not isinstance(profiles_raw, list) or not profiles_raw:
        raise ProfileValidationError(
            "population must declare a non-empty 'profiles' array"
        )

    # Determine the configured loci. Explicit 'loci' wins; otherwise infer from
    # the first profile and require every profile to match it exactly.
    if "loci" in data:
        if not isinstance(data["loci"], list) or not data["loci"]:
            raise ProfileValidationError("'loci' must be a non-empty array")
        loci = tuple(_normalize_locus(x) for x in data["loci"])
        if len(set(loci)) != len(loci):
            raise ProfileValidationError(f"'loci' contains duplicates: {loci}")
    else:
        first = profiles_raw[0]
        if not isinstance(first, dict) or not isinstance(first.get("resolution"), dict):
            raise ProfileValidationError(
                "cannot infer loci: first profile lacks a 'resolution' table"
            )
        loci = tuple(_normalize_locus(k) for k in first["resolution"])

    profiles = tuple(build_profile(p, loci=loci) for p in profiles_raw)

    ids = [p.id for p in profiles]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        raise ProfileValidationError(f"duplicate profile ids: {dupes}")

    total = math.fsum(p.weight for p in profiles)
    if total <= 0:
        raise ProfileValidationError(
            "at least one profile must have positive weight (total weight is 0)"
        )
    if abs(total - 1.0) > _WEIGHT_SUM_TOLERANCE:
        raise ProfileValidationError(
            f"profile weights must sum to 1 (within {_WEIGHT_SUM_TOLERANCE}); "
            f"got {total}"
        )

    return PopulationTypingModel(name=name, loci=loci, profiles=profiles)
