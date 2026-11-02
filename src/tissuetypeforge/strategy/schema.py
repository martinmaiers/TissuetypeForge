"""Machine-readable typing-strategy schema (spec §3).

A strategy is a configuration over a common observation framework — not a
hard-coded conversion function (spec §3). It specifies, per locus, whether the
locus is typed, the assay class, the target resolution, the reference database
version, and (in later milestones) sequence regions, ambiguity rules, failure
and error probabilities.

Milestone 1 implements the *idealized* resolution modes only. Assay classes and
resolution modes that require reference sequence or serology data are declared
in the enums so strategies can be written against the final vocabulary, but the
observation framework raises an explicit ``NotImplementedError`` for them in M1
rather than faking a transformation (spec §12).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class Resolution(StrEnum):
    """Target resolution of a typing observation."""

    #: Report the true allele verbatim (idealized, M1).
    ALLELE = "allele"
    #: Truncate to the first nomenclature field (idealized, M1).
    FIRST_FIELD = "first_field"
    #: Truncate to the first two nomenclature fields (idealized, M1).
    SECOND_FIELD = "second_field"
    #: Serologic specificity (requires versioned mapping — M2).
    SEROLOGY = "serology"
    #: Intermediate resolution from interrogated regions (M2).
    INTERMEDIATE = "intermediate"
    #: Resolve to the extent complete reference sequence permits (M2).
    FULL_GENE = "full_gene"

    @property
    def is_idealized_m1(self) -> bool:
        return self in {Resolution.ALLELE, Resolution.FIRST_FIELD, Resolution.SECOND_FIELD}


class AssayClass(StrEnum):
    """Broad class of typing assay (spec §3)."""

    #: Idealized nomenclature truncation; no assay physics modeled (M1).
    IDEALIZED = "idealized"
    #: Serologic typing (M2).
    SEROLOGY = "serology"
    #: Sequence-specific oligonucleotide / probe-based (M2).
    SSO = "sso"
    #: Sequence-based typing / NGS over interrogated regions (M2).
    SBT = "sbt"


class StrategyValidationError(ValueError):
    """Raised when a strategy configuration is invalid."""


@dataclass(frozen=True)
class LocusConfig:
    """Per-locus observation configuration.

    M1 fields plus M2 additions (Decision 4):

    * ``reference_name`` selects a versioned reference artifact from the registry
      (sequence column model or serology model) for sequence/serology assays.
    * ``interrogated_features`` names the features an assay reads, e.g.
      ``("exon2", "exon3")`` — used by the sequence engine to pick columns.
    * ``reference_policy`` is the per-locus strict/permissive/sequence_supplied
      handling of incomplete reference (overrides the strategy default).
    * ``serology_level`` is ``broad`` or ``split`` for serology assays.
    * ``observed_epitopes`` is the explicitly declared list of public epitopes
      (e.g. ``("Bw4", "Bw6")``) the assay interrogates (Decision 1). Only these
      may participate in a serologic observation.

    ``assay_params`` holds optional assay-specific parameters (spec §3).
    """

    locus: str
    typed: bool = True
    assay_class: AssayClass = AssayClass.IDEALIZED
    resolution: Resolution = Resolution.SECOND_FIELD
    reference_database: str | None = None
    reference_name: str | None = None
    interrogated_features: tuple[str, ...] = ()
    reference_policy: str | None = None
    serology_level: str = "split"
    observed_epitopes: tuple[str, ...] = ()
    sequence_regions: tuple[str, ...] = ()
    ambiguity_rules: str | None = None
    failure_probability: float = 0.0
    error_probability: float = 0.0
    assay_params: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.locus:
            raise StrategyValidationError("LocusConfig.locus must be non-empty")
        for name, p in (
            ("failure_probability", self.failure_probability),
            ("error_probability", self.error_probability),
        ):
            if not (0.0 <= p <= 1.0):
                raise StrategyValidationError(
                    f"{name} for locus {self.locus} must be in [0, 1], got {p}"
                )
        if self.serology_level not in ("broad", "split"):
            raise StrategyValidationError(
                f"serology_level for {self.locus} must be 'broad' or 'split', "
                f"got {self.serology_level!r}"
            )


@dataclass(frozen=True)
class PhaseGroup:
    """A set of loci whose interlocus phase the assay retains (spec §6).

    Loci not sharing any phase group are reported with no interlocus phase.
    A single-locus phase group is meaningless and rejected.
    """

    loci: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.loci) < 2:
            raise StrategyValidationError(
                f"A phase group must contain at least two loci, got {self.loci}"
            )
        if len(set(self.loci)) != len(self.loci):
            raise StrategyValidationError(f"Phase group has duplicate loci: {self.loci}")

    @property
    def locus_set(self) -> frozenset[str]:
        return frozenset(self.loci)


@dataclass(frozen=True)
class Strategy:
    """A complete typing strategy configuration.

    ``loci`` maps locus key (e.g. ``HLA-A`` or ``A``) to its config.
    ``phase_groups`` declares interlocus phase retention. A locus with no
    explicit config falls back to ``default_locus`` when it appears in truth.
    """

    name: str
    loci: dict[str, LocusConfig]
    phase_groups: tuple[PhaseGroup, ...] = ()
    reference_database: str | None = None
    description: str = ""
    default_locus: LocusConfig | None = None
    #: Default reference policy for sequence/serology loci (M2): one of
    #: 'strict' (default), 'permissive', 'sequence_supplied'. A locus may
    #: override via LocusConfig.reference_policy.
    reference_policy: str = "strict"
    #: Optional resource limits (M2): mapping of limit name -> int|None, e.g.
    #: {"max_universe": 100000, "max_alternatives": 10000}. Empty = defaults.
    resource_limits: dict[str, int | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise StrategyValidationError("Strategy.name must be non-empty")
        if self.reference_policy not in ("strict", "permissive", "sequence_supplied"):
            raise StrategyValidationError(
                f"reference_policy must be strict/permissive/sequence_supplied, "
                f"got {self.reference_policy!r}"
            )
        # A locus may appear in at most one phase group (disjoint groups).
        seen: set[str] = set()
        for pg in self.phase_groups:
            overlap = seen & pg.locus_set
            if overlap:
                raise StrategyValidationError(
                    f"Loci {overlap} appear in more than one phase group"
                )
            seen |= pg.locus_set

    def config_for(self, locus: str) -> LocusConfig | None:
        """Return the config for ``locus``, falling back to the default.

        Matches both the full ``HLA-A`` form and a bare ``A`` key so strategies
        can be written compactly. Returns ``None`` only if no config and no
        default exist (locus is then treated as not typed).
        """
        if locus in self.loci:
            return self.loci[locus]
        short = locus.split("-", 1)[-1]
        if short in self.loci:
            return self.loci[short]
        if self.default_locus is not None:
            # materialize a default bound to this locus name
            d = self.default_locus
            return LocusConfig(
                locus=locus,
                typed=d.typed,
                assay_class=d.assay_class,
                resolution=d.resolution,
                reference_database=d.reference_database,
                reference_name=d.reference_name,
                interrogated_features=d.interrogated_features,
                reference_policy=d.reference_policy,
                serology_level=d.serology_level,
                observed_epitopes=d.observed_epitopes,
                sequence_regions=d.sequence_regions,
                ambiguity_rules=d.ambiguity_rules,
                failure_probability=d.failure_probability,
                error_probability=d.error_probability,
                assay_params=dict(d.assay_params),
            )
        return None

    def phase_group_for(self, locus: str) -> PhaseGroup | None:
        short = locus.split("-", 1)[-1]
        for pg in self.phase_groups:
            if locus in pg.loci or short in pg.loci:
                return pg
        return None
