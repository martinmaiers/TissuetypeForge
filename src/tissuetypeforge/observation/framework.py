"""Observation framework core (spec §2, §5, §6; Milestone 1).

Pipeline:

    truth GL String
      -> parse -> Diplotype
      -> decompose to per-locus observed genotypes (coverage + resolution)
      -> regroup by phase-retention rules -> observed Diplotype
      -> serialize -> observed GL String + ObservationResult

Phase handling (spec §6) is the subtle part:

* Phase is only *retained*, never *invented*. Two loci may be reported in one
  phased block only if (a) truth phased them together in one GenotypeBlock and
  (b) the strategy declares a phase group containing both. Otherwise each locus
  is reported as an independent block, asserting no interlocus phase.
* Breaking phase therefore cannot create new phase information: it only ever
  splits a block into per-locus blocks (QA §10.2).

In M1 there is no stochastic error, so ``observed == expected_observation`` and
the error category is ``NONE``. Assay ambiguity / information loss is recorded
in ``loss_record`` and is never labeled mistyping (spec §7).
"""

from __future__ import annotations

from dataclasses import dataclass

from tissuetypeforge.glstring import (
    AlleleList,
    Genotype,
    GenotypeList,
    LocusBlockList,
    PhasedGenotype,
    parse_blocks,
    serialize,
)
from tissuetypeforge.observation.joint_phase import (
    enumerate_joint,
    joint_to_genotype_list,
)
from tissuetypeforge.observation.policy import (
    DEFAULT_LIMITS,
    ReferencePolicy,
    ResourceLimits,
)
from tissuetypeforge.observation.resolution import resolve_allele_list
from tissuetypeforge.observation.sequence_engine import enumerate_locus, to_genotype_list
from tissuetypeforge.reference.model import FeatureKey
from tissuetypeforge.reference.registry import ReferenceRegistry
from tissuetypeforge.result import (
    ErrorCategory,
    LocusObservation,
    ObservationMetadata,
    ObservationResult,
)
from tissuetypeforge.serology.backend import observe_serology
from tissuetypeforge.serology.model import Epitope, ResolutionLevel
from tissuetypeforge.strategy.schema import AssayClass, LocusConfig, Strategy


class ObservationError(ValueError):
    """Raised when truth cannot be observed under the given strategy."""


def _parse_feature(token: str) -> FeatureKey:
    """Parse a feature token like ``exon2`` / ``intron1`` into a FeatureKey."""
    kind_chars = token.rstrip("0123456789")
    idx = token[len(kind_chars):]
    if not kind_chars or not idx:
        raise ObservationError(f"malformed interrogated feature token {token!r}")
    return FeatureKey(kind_chars, int(idx))


def _features(config: LocusConfig) -> frozenset[FeatureKey]:
    return frozenset(_parse_feature(t) for t in config.interrogated_features)


def _resource_limits(strategy: Strategy) -> ResourceLimits:
    rl = strategy.resource_limits
    if not rl:
        return DEFAULT_LIMITS
    base = ResourceLimits()
    return ResourceLimits(
        max_universe=rl.get("max_universe", base.max_universe),
        max_candidate_pairs=rl.get("max_candidate_pairs", base.max_candidate_pairs),
        max_alternatives=rl.get("max_alternatives", base.max_alternatives),
        max_joint_molecules=rl.get("max_joint_molecules", base.max_joint_molecules),
    )


def _reference_policy(strategy: Strategy, config: LocusConfig) -> ReferencePolicy:
    name = config.reference_policy or strategy.reference_policy
    return ReferencePolicy(name)


@dataclass
class _LocusTruth:
    """The true genotype at a single locus, plus its truth-phase identity.

    ``copies`` has one AlleleList per gene copy (two for a diploid locus).
    ``truth_block_id`` identifies the original truth locus block this locus came
    from, so phase can only be retained among loci truly in phase.
    """

    locus: str
    copies: list[AlleleList]
    truth_block_id: int


def _decompose_truth(blocks: LocusBlockList) -> tuple[list[_LocusTruth], list[str]]:
    """Flatten truth into per-locus genotypes, preserving phase identity.

    Rejects genotype-ambiguous truth in M1: truth is assumed known and
    unambiguous (the simulator *produces* ambiguity, it does not consume it).
    Returns the per-locus truths and the locus order as they appear in truth.
    """
    locus_truths: list[_LocusTruth] = []
    locus_order: list[str] = []
    for block_id, genlist in enumerate(blocks.blocks):
        if genlist.is_ambiguous:
            raise ObservationError(
                "Genotype-ambiguous truth is not supported; provide known truth "
                "as an unambiguous GL String"
            )
        genotype = genlist.genotypes[0]
        # Ordered loci across the (possibly phased) block: the loci of copy 0.
        loci = genotype.copies[0].loci
        # For each locus position within the (possibly phased) block, gather the
        # AlleleList from each copy at that position.
        for pos, locus in enumerate(loci):
            copies_at_locus = [copy.segments[pos] for copy in genotype.copies]
            locus_truths.append(
                _LocusTruth(
                    locus=locus,
                    copies=copies_at_locus,
                    truth_block_id=block_id,
                )
            )
            locus_order.append(locus)
    return locus_truths, locus_order


def _observe_locus(
    truth: _LocusTruth, config: LocusConfig
) -> tuple[list[AlleleList], list[str]]:
    """Apply idealized resolution to each copy at one locus.

    Returns the observed per-copy AlleleLists and any loss notes.
    """
    notes: list[str] = []
    observed_copies: list[AlleleList] = []
    for seg in truth.copies:
        resolved = resolve_allele_list(seg, config.resolution)
        observed_copies.append(resolved)
    # Record information loss: did truncation drop fields?
    for orig, obs in zip(truth.copies, observed_copies, strict=True):
        o_names = {a.name for a in orig.alleles}
        n_names = {a.name for a in obs.alleles}
        if o_names != n_names:
            notes.append(
                f"{truth.locus}: resolution '{config.resolution.value}' reduced "
                f"{sorted(o_names)} -> {sorted(n_names)}"
            )
    return observed_copies, notes


def observe(
    truth_gl: str,
    strategy: Strategy,
    *,
    seed: int | None = None,
    registry: ReferenceRegistry | None = None,
) -> ObservationResult:
    """Observe ``truth_gl`` under ``strategy`` and return a structured result.

    Dispatches per locus by ``assay_class`` (Decision 4):

    * ``idealized`` -> M1 nomenclature-truncation path (unchanged, no reference);
    * ``sbt`` / ``sso`` -> sequence engine over a registry-selected column model;
    * ``serology`` -> serology backend (antigen report + compatible-DNA GL).

    Sequence and serology loci require a ``registry`` providing the named
    reference artifacts. ``seed`` is recorded but unused (no stochastic error in
    M2). Backward-compatible with M1: an all-idealized strategy needs no
    ``registry`` and produces identical output.
    """
    truth_blocks = parse_blocks(truth_gl)
    locus_truths, _ = _decompose_truth(truth_blocks)
    limits = _resource_limits(strategy)

    typed_loci: list[str] = []
    missing_loci: list[str] = []
    loss_record: list[str] = []
    locus_observations: list[LocusObservation] = []
    reference_dbs: set[str] = set()
    antigen_reports: list[object] = []

    # idealized loci go through M1 phase-regrouping; M2 loci produce standalone
    # per-locus blocks appended directly (independent ^ groups).
    observed_by_locus: dict[str, list[AlleleList]] = {}
    truth_block_by_locus: dict[str, int] = {}
    m2_blocks: list[GenotypeList] = []

    lt_by_locus = {lt.locus: lt for lt in locus_truths}
    retained_joint: list[tuple[str, ...]] = []
    consumed = _observe_joint_groups(
        strategy, lt_by_locus, registry, limits,
        m2_blocks, typed_loci, loss_record, locus_observations,
        reference_dbs, retained_joint,
    )

    for lt in locus_truths:
        if lt.locus in consumed:
            continue
        config = strategy.config_for(lt.locus)
        if config is None or not config.typed:
            missing_loci.append(lt.locus)
            locus_observations.append(
                LocusObservation(
                    locus=lt.locus, typed=False, resolution=None,
                    assay_class=None, reference_database=None,
                    notes=("locus not typed by strategy (missing result)",),
                )
            )
            continue

        if config.assay_class is AssayClass.IDEALIZED:
            observed_copies, notes = _observe_locus(lt, config)
            observed_by_locus[lt.locus] = observed_copies
            truth_block_by_locus[lt.locus] = lt.truth_block_id
            loss_record.extend(notes)
            db = config.reference_database
        else:
            block, notes, db, report = _observe_m2_locus(
                lt, config, strategy, registry, limits
            )
            m2_blocks.append(block)
            loss_record.extend(notes)
            if report is not None:
                antigen_reports.append(report)

        typed_loci.append(lt.locus)
        if db:
            reference_dbs.add(db)
        locus_observations.append(
            LocusObservation(
                locus=lt.locus, typed=True,
                resolution=config.resolution.value,
                assay_class=config.assay_class.value,
                reference_database=db,
                notes=tuple(notes),
            )
        )

    idealized_groups, retained_phase = _assemble_groups(
        [locus for locus in typed_loci if locus in observed_by_locus],
        observed_by_locus, truth_block_by_locus, strategy,
    )
    observed_groups = list(idealized_groups) + m2_blocks
    all_retained = list(retained_phase) + retained_joint

    if not observed_groups:
        observed_gl = ""
    else:
        observed_blocks = LocusBlockList(tuple(observed_groups))
        observed_gl = serialize(observed_blocks)

    metadata = ObservationMetadata(
        strategy_name=strategy.name,
        reference_databases=tuple(sorted(reference_dbs)),
        typed_loci=tuple(typed_loci),
        missing_loci=tuple(missing_loci),
        retained_phase_groups=tuple(all_retained),
        locus_observations=tuple(locus_observations),
        seed=seed,
    )

    return ObservationResult(
        true_diplotype=truth_gl,
        expected_observation=observed_gl,
        observed=observed_gl,  # M2: no stochastic error
        error_category=ErrorCategory.NONE,
        metadata=metadata,
        error_events=(),
        loss_record=tuple(loss_record),
        antigen_reports=tuple(antigen_reports),
    )


def _observe_m2_locus(lt, config, strategy, registry, limits):
    """Dispatch a single sequence/serology locus; return (block, notes, db, report)."""
    if registry is None:
        raise ObservationError(
            f"{config.assay_class.value} assay at {lt.locus} needs a reference "
            f"registry; call observe(..., registry=...)"
        )
    true_pair = _true_pair(lt)

    if config.assay_class in (AssayClass.SBT, AssayClass.SSO):
        if not config.reference_name:
            raise ObservationError(
                f"sequence assay at {lt.locus} requires reference_name"
            )
        model = registry.column_model(config.reference_name, lt.locus)
        cand = enumerate_locus(
            model, true_pair, _features(config),
            reference_policy=_reference_policy(strategy, config),
            limits=limits,
        )
        block = to_genotype_list(cand, limits=limits)
        notes = [
            f"{lt.locus}: {config.assay_class.value} over "
            f"{sorted(config.interrogated_features)} -> "
            f"{len(cand.all_pairs)} candidate genotype(s)"
        ]
        if cand.indeterminate:
            notes.append(
                f"{lt.locus}: {len(cand.indeterminate)} indeterminate candidate(s) "
                f"from incomplete reference (permissive)"
            )
        return block, notes, model.database_version, None

    if config.assay_class is AssayClass.SEROLOGY:
        if not config.reference_name:
            raise ObservationError(
                f"serology assay at {lt.locus} requires reference_name"
            )
        model = registry.serology_model(config.reference_name, lt.locus)
        level = ResolutionLevel(config.serology_level)
        obs = observe_serology(
            model, true_pair, level,
            observed_epitopes=frozenset(Epitope(e) for e in config.observed_epitopes),
            limits=limits,
        )
        notes = [
            f"{lt.locus}: serology {level.value} -> "
            f"antigens {sorted(s.name for s in obs.report.specificities)}; "
            f"{len(obs.compatible_pairs)} compatible genotype(s)"
        ]
        return obs.compatible_genotypes, notes, model.database_version, obs.report

    raise ObservationError(f"unsupported assay_class {config.assay_class} at {lt.locus}")


def _observe_joint_groups(
    strategy, lt_by_locus, registry, limits,
    m2_blocks, typed_loci, loss_record, locus_observations,
    reference_dbs, retained_joint,
) -> set[str]:
    """Pre-pass: enumerate declared phase groups whose loci are sequence-typed as
    whole-molecule joint observations (D8, main §6.4.1). Returns consumed loci.

    A phase group participates only if (a) all its loci are present in truth,
    typed, and sequence-class (sbt/sso), and (b) all its loci share one truth
    block (truly phased). Otherwise the group is left to the per-locus path
    (phase is only retained where observed - Decision 3 / §6.4 Rule A).
    """
    consumed: set[str] = set()
    for pg in strategy.phase_groups:
        loci: list[str] = []
        configs: list[LocusConfig] = []
        ok = True
        for raw in pg.loci:
            match = None
            for name in lt_by_locus:
                short = name.split("-", 1)[-1]
                if name == raw or short == raw:
                    match = name
                    break
            if match is None:
                ok = False
                break
            cfg = strategy.config_for(match)
            if cfg is None or not cfg.typed or cfg.assay_class not in (
                AssayClass.SBT, AssayClass.SSO
            ):
                ok = False
                break
            loci.append(match)
            configs.append(cfg)
        if not ok or len(loci) < 2:
            continue
        block_ids = {lt_by_locus[x].truth_block_id for x in loci}
        if len(block_ids) != 1:
            # truth did not phase them -> do not invent phase; per-locus path
            continue
        if registry is None:
            raise ObservationError(f"joint-phase group {loci} needs a registry")

        loci_t = tuple(loci)
        models = {
            x: registry.column_model(configs[i].reference_name, x)
            for i, x in enumerate(loci)
        }
        features = {x: _features(configs[i]) for i, x in enumerate(loci)}
        hap0 = tuple(lt_by_locus[x].copies[0].alleles[0].name for x in loci)
        hap1 = tuple(lt_by_locus[x].copies[1].alleles[0].name for x in loci)
        cand = enumerate_joint(
            models, loci_t, (hap0, hap1), features,
            reference_policy=ReferencePolicy(strategy.reference_policy),
            limits=limits,
        )
        m2_blocks.append(joint_to_genotype_list(cand, limits=limits))
        retained_joint.append(loci_t)
        for i, x in enumerate(loci):
            consumed.add(x)
            typed_loci.append(x)
            reference_dbs.add(models[x].database_version)
            locus_observations.append(
                LocusObservation(
                    locus=x, typed=True,
                    resolution=configs[i].resolution.value,
                    assay_class=configs[i].assay_class.value,
                    reference_database=models[x].database_version,
                    notes=(f"{x}: joint phase observation over {loci_t}; "
                           f"{len(cand.confirmed)} joint candidate(s)",),
                )
            )
        loss_record.append(
            f"joint phase group {loci_t}: {len(cand.confirmed)} candidate diplotype(s)"
        )
    return consumed


def _true_pair(lt) -> tuple[str, str]:
    """Extract the two true allele names at a single-locus truth."""
    names = []
    for seg in lt.copies:
        if len(seg.alleles) != 1:
            raise ObservationError(
                f"{lt.locus}: truth must be an unambiguous allele per copy"
            )
        names.append(seg.alleles[0].name)
    if len(names) != 2:
        raise ObservationError(f"{lt.locus}: expected a diploid pair, got {names}")
    return (names[0], names[1]) if names[0] <= names[1] else (names[1], names[0])


def _assemble_groups(
    typed_loci: list[str],
    observed_by_locus: dict[str, list[AlleleList]],
    truth_block_by_locus: dict[str, int],
    strategy: Strategy,
) -> tuple[list[GenotypeList], list[tuple[str, ...]]]:
    """Regroup typed loci into observed locus blocks per phase-retention rules.

    Two loci are reported in one phased block iff they were in the same truth
    block (truly in phase) AND the strategy declares a phase group containing
    both. Phase is never invented (spec §6). Returns the observed blocks in
    truth order and the list of phase groups actually retained.
    """
    # Partition typed loci into phased bundles.
    # Key: (truth_block_id, phase_group_key) where phase_group_key is the frozen
    # locus set of the declared phase group, or a unique marker if none.
    bundles: list[list[str]] = []
    bundle_index: dict[tuple[int, object], int] = {}
    retained_phase: list[tuple[str, ...]] = []

    for locus in typed_loci:
        block_id = truth_block_by_locus[locus]
        pg = strategy.phase_group_for(locus)
        if pg is None:
            # No declared phase group -> independent locus, unique bundle.
            key: tuple[int, object] = (block_id, ("__nophase__", locus))
        else:
            key = (block_id, pg.locus_set)
        if key not in bundle_index:
            bundle_index[key] = len(bundles)
            bundles.append([])
        bundles[bundle_index[key]].append(locus)

    groups: list[GenotypeList] = []
    for bundle in bundles:
        if len(bundle) > 1:
            retained_phase.append(tuple(bundle))
        n_copies = len(observed_by_locus[bundle[0]])
        # Build each copy as a ~-ordered run across the bundle's loci.
        copies: list[PhasedGenotype] = []
        for copy_idx in range(n_copies):
            segments = tuple(observed_by_locus[locus][copy_idx] for locus in bundle)
            copies.append(PhasedGenotype(segments))
        genotype = Genotype(tuple(copies))
        groups.append(GenotypeList((genotype,)))

    return groups, retained_phase
