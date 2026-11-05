"""TissueTypeForge: simulate HLA typing observations from known genetic truth.

Simplified observation model with population-level joint typing-profile
distributions.

Core concepts:
    ResolutionMode                 controlled vocabulary: 0-field/1-field/2-field
    TypingProfile                  a complete per-locus resolution assignment
    PopulationTypingModel          categorical distribution over profiles
    observe_genotype               observe one donor under one profile
    simulate_cohort                seeded cohort simulation + reproducibility manifest
    parse_gl_string / serialize    general-purpose GL String round-trip

Observations are *unphased*: no interlocus phase is asserted. A 0-field locus is
a missing result (omitted), never a fabricated homozygous call.
"""

from tissuetypeforge.config import load_population, load_population_file
from tissuetypeforge.glstring import parse_gl_string, serialize
from tissuetypeforge.observe import (
    DonorObservation,
    LocusObservation,
    ObservationError,
    observe_genotype,
)
from tissuetypeforge.population import (
    PopulationTypingModel,
    ProfileValidationError,
    TypingProfile,
    build_population,
    build_profile,
)
from tissuetypeforge.resolution import (
    ResolutionError,
    ResolutionMode,
    parse_resolution,
)
from tissuetypeforge.simulate import (
    CohortResult,
    ReproducibilityManifest,
    sample_profile,
    simulate_cohort,
)

__all__ = [
    "CohortResult",
    "DonorObservation",
    "LocusObservation",
    "ObservationError",
    "PopulationTypingModel",
    "ProfileValidationError",
    "ReproducibilityManifest",
    "ResolutionError",
    "ResolutionMode",
    "TypingProfile",
    "build_population",
    "build_profile",
    "load_population",
    "load_population_file",
    "observe_genotype",
    "parse_gl_string",
    "parse_resolution",
    "sample_profile",
    "serialize",
    "simulate_cohort",
]

__version__ = "0.3.0"
