"""Tests for profiles, population config, validation, sampling, reproducibility.

Covers REQ 7 items 4 (joint sampling), 5 (proportions converge), 6 (determinism),
7 (invalid configs rejected).
"""

from __future__ import annotations

import math

import pytest

from tissuetypeforge import (
    ProfileValidationError,
    build_population,
    load_population,
    simulate_cohort,
)

# Example population from the requirement (weights 0.2 / 0.5 / 0.3).
POP_DATA = {
    "name": "registry-x",
    "loci": ["A", "B", "C", "DRB1"],
    "profiles": [
        {
            "id": "legacy",
            "weight": 0.20,
            "resolution": {"A": "1-field", "B": "1-field", "C": "0-field", "DRB1": "1-field"},
        },
        {
            "id": "intermediate",
            "weight": 0.50,
            "resolution": {"A": "2-field", "B": "2-field", "C": "1-field", "DRB1": "2-field"},
        },
        {
            "id": "modern",
            "weight": 0.30,
            "resolution": {"A": "2-field", "B": "2-field", "C": "2-field", "DRB1": "2-field"},
        },
    ],
}

TRUTH = (
    "HLA-A*02:01:01:01+HLA-A*03:01:01:01^"
    "HLA-B*07:02:01+HLA-B*08:01:01^"
    "HLA-C*07:02:01+HLA-C*07:01:01^"
    "HLA-DRB1*15:01:01+HLA-DRB1*03:01:01"
)


# --- REQ 7.7: invalid configurations rejected -----------------------------


def test_duplicate_profile_ids_rejected() -> None:
    data = {
        "loci": ["A"],
        "profiles": [
            {"id": "x", "weight": 0.5, "resolution": {"A": "1-field"}},
            {"id": "x", "weight": 0.5, "resolution": {"A": "2-field"}},
        ],
    }
    with pytest.raises(ProfileValidationError, match="duplicate profile ids"):
        build_population(data)


def test_negative_weight_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": -0.1, "resolution": {"A": "1-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="nonnegative"):
        build_population(data)


def test_nonfinite_weight_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": float("inf"), "resolution": {"A": "1-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="finite"):
        build_population(data)


def test_zero_total_weight_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": 0.0, "resolution": {"A": "1-field"}},
        {"id": "y", "weight": 0.0, "resolution": {"A": "2-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="positive weight"):
        build_population(data)


def test_weights_not_summing_to_one_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": 0.2, "resolution": {"A": "1-field"}},
        {"id": "y", "weight": 0.2, "resolution": {"A": "2-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="sum to 1"):
        build_population(data)


def test_profile_missing_a_configured_locus_rejected() -> None:
    data = {"loci": ["A", "B"], "profiles": [
        {"id": "x", "weight": 1.0, "resolution": {"A": "1-field"}},  # missing B
    ]}
    with pytest.raises(ProfileValidationError, match="missing configured loci"):
        build_population(data)


def test_unknown_locus_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": 1.0, "resolution": {"A": "1-field", "Z": "1-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="unknown loci"):
        build_population(data)


def test_unsupported_resolution_mode_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"id": "x", "weight": 1.0, "resolution": {"A": "3-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="Unsupported resolution"):
        build_population(data)


def test_empty_profiles_rejected() -> None:
    with pytest.raises(ProfileValidationError, match="non-empty 'profiles'"):
        build_population({"loci": ["A"], "profiles": []})


def test_missing_id_rejected() -> None:
    data = {"loci": ["A"], "profiles": [
        {"weight": 1.0, "resolution": {"A": "1-field"}},
    ]}
    with pytest.raises(ProfileValidationError, match="non-empty 'id'"):
        build_population(data)


def test_valid_population_accepted() -> None:
    pop = build_population(POP_DATA)
    assert pop.name == "registry-x"
    assert pop.loci == ("A", "B", "C", "DRB1")
    assert [p.id for p in pop.profiles] == ["legacy", "intermediate", "modern"]
    assert math.isclose(sum(pop.normalized_weights()), 1.0)


def test_loci_inferred_from_first_profile_when_omitted() -> None:
    data = {"profiles": [
        {"id": "x", "weight": 1.0, "resolution": {"A": "1-field", "B": "2-field"}},
    ]}
    pop = build_population(data)
    assert set(pop.loci) == {"A", "B"}


def test_load_population_accepts_hla_prefixed_loci() -> None:
    data = {"loci": ["HLA-A"], "profiles": [
        {"id": "x", "weight": 1.0, "resolution": {"HLA-A": "1-field"}},
    ]}
    pop = load_population(data)
    assert pop.loci == ("A",)


# --- REQ 7.4: joint (not independent per-locus) sampling ------------------


def test_sampling_is_joint_not_independent_per_locus() -> None:
    # legacy: C=0-field; modern: C=2-field, and legacy has A=1-field while
    # modern has A=2-field. Under JOINT sampling, a donor whose C is missing
    # (legacy) must ALSO show A at 1-field — the per-locus choices are perfectly
    # correlated through the profile. Independent per-locus sampling would mix
    # them. We assert the correlation holds for every donor.
    pop = build_population(POP_DATA)
    result = simulate_cohort([TRUTH] * 300, pop, seed=7)
    for obs in result.observations:
        a_one_field = "HLA-A*02+" in obs.observed or obs.observed.startswith("HLA-A*02^")
        c_missing = "HLA-C" in obs.missing_loci
        if obs.profile_id == "legacy":
            # legacy => A at 1-field AND C missing, jointly.
            assert a_one_field and c_missing
        else:
            # non-legacy => A at 2-field (never 1-field) and C present.
            assert not c_missing


def test_assignment_distribution_matches_profile_ids() -> None:
    pop = build_population(POP_DATA)
    result = simulate_cohort([TRUTH] * 50, pop, seed=1)
    # every assignment is one of the configured profile ids
    assert set(result.manifest.profile_assignments) <= {"legacy", "intermediate", "modern"}
    assert result.manifest.n_donors == 50


# --- REQ 7.5: proportions approach configured probabilities ---------------


def test_observed_proportions_converge_for_large_cohort() -> None:
    pop = build_population(POP_DATA)
    n = 20000
    result = simulate_cohort([TRUTH] * n, pop, seed=12345)
    counts = result.manifest.assigned_profile_counts
    expected = {"legacy": 0.20, "intermediate": 0.50, "modern": 0.30}
    for pid, exp in expected.items():
        observed = counts[pid] / n
        assert abs(observed - exp) < 0.02, f"{pid}: {observed} vs {exp}"


# --- REQ 7.6: deterministic reproducibility -------------------------------


def test_same_seed_same_outputs() -> None:
    pop = build_population(POP_DATA)
    r1 = simulate_cohort([TRUTH] * 500, pop, seed=42)
    r2 = simulate_cohort([TRUTH] * 500, pop, seed=42)
    assert r1.manifest.profile_assignments == r2.manifest.profile_assignments
    assert [o.observed for o in r1.observations] == [o.observed for o in r2.observations]


def test_different_seed_differs() -> None:
    pop = build_population(POP_DATA)
    r1 = simulate_cohort([TRUTH] * 500, pop, seed=1)
    r2 = simulate_cohort([TRUTH] * 500, pop, seed=2)
    assert r1.manifest.profile_assignments != r2.manifest.profile_assignments


def test_manifest_records_profile_ids_and_seed() -> None:
    pop = build_population(POP_DATA)
    result = simulate_cohort([TRUTH] * 10, pop, seed=99, donor_ids=[f"D{i}" for i in range(10)])
    m = result.manifest
    assert m.seed == 99
    assert m.population_name == "registry-x"
    assert set(m.profile_weights) == {"legacy", "intermediate", "modern"}
    assert math.isclose(sum(m.profile_weights.values()), 1.0)
    assert sum(m.assigned_profile_counts.values()) == 10
    assert len(m.profile_assignments) == 10
    assert result.donor_ids == tuple(f"D{i}" for i in range(10))


def test_donor_ids_length_mismatch_rejected() -> None:
    pop = build_population(POP_DATA)
    with pytest.raises(ValueError, match="donor_ids length"):
        simulate_cohort([TRUTH, TRUTH], pop, seed=1, donor_ids=["only-one"])
