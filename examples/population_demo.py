"""Population typing-profile simulation demo (TissueTypeForge).

Run with:  uv run python examples/population_demo.py

Demonstrates the simplified observation model:
  * a controlled resolution vocabulary (0-field / 1-field / 2-field);
  * a population as a JOINT distribution over complete typing profiles;
  * seeded, reproducible cohort simulation;
  * a reproducibility manifest recording profile assignments.
"""

from __future__ import annotations

from pathlib import Path

from tissuetypeforge import load_population_file, observe_genotype, simulate_cohort

POP_FILE = Path(__file__).resolve().parents[1] / "populations" / "registry_example.toml"

# Known true genotypes (e.g. from DiplotypeForge or a real source). The genetic
# composition is supplied independently of the typing model.
TRUTHS = [
    "HLA-A*02:01:01:01+HLA-A*03:01:01:01^"
    "HLA-C*07:02:01+HLA-C*07:01:01^"
    "HLA-B*07:02:01+HLA-B*08:01:01^"
    "HLA-DRB1*15:01:01+HLA-DRB1*03:01:01",
    "HLA-A*01:01:01:01+HLA-A*02:01:01:01^"
    "HLA-C*03:04:01+HLA-C*06:02:01^"
    "HLA-B*35:01:01+HLA-B*44:02:01^"
    "HLA-DRB1*04:01:01+HLA-DRB1*07:01:01",
]


def _rule() -> None:
    print("-" * 72)


def main() -> None:
    pop = load_population_file(POP_FILE)
    print(f"Population: {pop.name}")
    print(f"Configured loci: {', '.join(pop.loci)}")
    print("Profiles (joint resolution across loci):")
    for profile, w in zip(pop.profiles, pop.normalized_weights(), strict=True):
        modes = ", ".join(f"{loc}={profile.resolution[loc].value}" for loc in pop.loci)
        print(f"  {profile.id:<12} p={w:.2f}  [{modes}]")
    _rule()

    # Apply each profile to one donor's truth to show the transformations.
    print("One donor's truth observed under each profile:")
    print(f"  truth: {TRUTHS[0]}")
    for profile in pop.profiles:
        obs = observe_genotype(TRUTHS[0], profile)
        missing = ", ".join(obs.missing_loci) or "(none)"
        print(f"  [{profile.id}] {obs.observed}")
        print(f"      missing: {missing}")
    _rule()

    # Seeded cohort simulation over a larger synthetic cohort.
    cohort = TRUTHS * 500  # 1000 donors
    donor_ids = [f"D{i:04d}" for i in range(len(cohort))]
    result = simulate_cohort(cohort, pop, seed=42, donor_ids=donor_ids)
    m = result.manifest
    print(f"Cohort simulation (seed={m.seed}, n={m.n_donors}):")
    for pid, count in m.assigned_profile_counts.items():
        print(f"  {pid:<12} assigned {count:5d}  ({count / m.n_donors:.3f}  "
              f"configured {m.profile_weights[pid]:.3f})")
    _rule()

    # Reproducibility: same seed -> identical assignments.
    again = simulate_cohort(cohort, pop, seed=42, donor_ids=donor_ids)
    identical = again.manifest.profile_assignments == m.profile_assignments
    print(f"Re-running with the same seed reproduces assignments exactly: {identical}")

    # No interlocus phase is ever asserted in observations.
    any_phase = any("~" in o.observed for o in result.observations)
    print(f"Any interlocus phase ('~') in observations: {any_phase}")


if __name__ == "__main__":
    main()
