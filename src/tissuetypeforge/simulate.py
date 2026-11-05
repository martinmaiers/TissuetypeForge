"""Seeded population/cohort simulation (REQ 2, 6).

For each donor the simulator:

1. samples exactly one typing profile from the population's categorical
   distribution (joint across loci — the whole profile is drawn at once, so
   resolution choices are correlated, not independent per locus);
2. applies that profile to the donor's true genotype;
3. produces the observed typing; and
4. records the assigned profile id.

Reproducibility (REQ 6): sampling uses a seeded :class:`random.Random`. Identical
inputs, configuration, and seed produce identical outputs. A
:class:`ReproducibilityManifest` captures the seed, population identity, profile
weights, and the per-donor profile assignments so a run can be audited and
replayed.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass, field

from tissuetypeforge.observe import DonorObservation, observe_genotype
from tissuetypeforge.population import PopulationTypingModel, TypingProfile


def sample_profile(
    population: PopulationTypingModel, rng: random.Random
) -> TypingProfile:
    """Draw exactly one profile from the population's categorical distribution.

    Uses ``random.Random.choices`` with the (unnormalized) profile weights — a
    single joint draw over complete profiles, never an independent per-locus
    draw (REQ 2).
    """
    weights = [p.weight for p in population.profiles]
    (chosen,) = rng.choices(population.profiles, weights=weights, k=1)
    return chosen


@dataclass(frozen=True)
class ReproducibilityManifest:
    """Audit record making a cohort run replayable (REQ 6)."""

    population_name: str
    seed: int | None
    loci: tuple[str, ...]
    #: profile id -> configured categorical probability (normalized).
    profile_weights: dict[str, float]
    #: profile id -> number of donors assigned that profile in this run.
    assigned_profile_counts: dict[str, int]
    #: ordered per-donor profile ids, in donor order (full assignment trace).
    profile_assignments: tuple[str, ...]
    n_donors: int


@dataclass(frozen=True)
class CohortResult:
    """Result of simulating observations for a cohort of donors."""

    observations: tuple[DonorObservation, ...]
    manifest: ReproducibilityManifest
    #: optional per-donor ids supplied by the caller, in donor order.
    donor_ids: tuple[str, ...] = field(default_factory=tuple)


def simulate_cohort(
    truth_genotypes: list[str],
    population: PopulationTypingModel,
    *,
    seed: int | None = None,
    donor_ids: list[str] | None = None,
) -> CohortResult:
    """Simulate observed typings for a cohort of donors (REQ 2, 6).

    ``truth_genotypes`` is one true genotype GL String per donor (the
    population's genetic composition — supplied independently of the typing
    model, preserving the REQ 3 separation of concerns). Each donor draws one
    profile and is observed under it.

    ``seed`` makes the run deterministic. ``donor_ids`` optionally labels donors;
    if given it must match ``truth_genotypes`` in length.
    """
    if donor_ids is not None and len(donor_ids) != len(truth_genotypes):
        raise ValueError(
            f"donor_ids length {len(donor_ids)} != truth_genotypes length "
            f"{len(truth_genotypes)}"
        )

    rng = random.Random(seed)
    observations: list[DonorObservation] = []
    assignments: list[str] = []

    for truth in truth_genotypes:
        profile = sample_profile(population, rng)
        assignments.append(profile.id)
        observations.append(observe_genotype(truth, profile))

    counts = Counter(assignments)
    weights = population.normalized_weights()
    manifest = ReproducibilityManifest(
        population_name=population.name,
        seed=seed,
        loci=population.loci,
        profile_weights={
            p.id: w for p, w in zip(population.profiles, weights, strict=True)
        },
        assigned_profile_counts={p.id: counts.get(p.id, 0) for p in population.profiles},
        profile_assignments=tuple(assignments),
        n_donors=len(truth_genotypes),
    )

    return CohortResult(
        observations=tuple(observations),
        manifest=manifest,
        donor_ids=tuple(donor_ids) if donor_ids is not None else (),
    )
