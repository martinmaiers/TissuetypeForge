"""Configuration model for the glnorm normalization harness.

A :class:`NormalizeConfig` captures every user-tunable knob so the pipeline is
fully reproducible: the same config applied to the same input yields
byte-identical output (see :mod:`glnorm.pipeline`).
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Valid py-ard reduction types (per the py-ard API). ``lgx`` is the default:
# a 2-field ARD reduction that is NOT naive field truncation.
VALID_REDUX_TYPES = ("G", "P", "lg", "lgx", "W", "exon", "U2", "S")

# The five classical loci, in the conventional serologic reporting order.
# These are matched against the ``PossibleLoci.loci`` names produced by
# glstring, which are prefixed (e.g. ``HLA-A``).
DEFAULT_LOCI = ("A", "C", "B", "DRB1", "DQB1")

DEFAULT_GL_PREFIX = "HLA-"


@dataclass(frozen=True)
class NormalizeConfig:
    """Immutable configuration for a normalization run.

    The three normalization operations are independent toggles:

    * ``normalize_nomenclature`` -- validate + canonicalize token form via
      glstring (never lossy, never expands ambiguity).
    * ``reduce_alleles`` -- OPTIONAL, lossy allele reduction via py-ard. Only
      applied when explicitly requested.
    * ``sort`` -- canonical ordering via gl_smartsort.
    """

    # --- Input parsing -----------------------------------------------------
    gl_column: str = "gl_string"
    donor_column: str | None = None
    input_format: str = "auto"  # auto | tsv | csv | lines

    # --- Locus configuration ----------------------------------------------
    loci: tuple[str, ...] = DEFAULT_LOCI
    locus_prefix: str = DEFAULT_GL_PREFIX
    #: If True, a record whose loci are not a subset of ``loci`` is an error.
    enforce_loci: bool = False

    # --- The three distinct operations (independently toggleable) ----------
    normalize_nomenclature: bool = True
    reduce_alleles: bool = False
    redux_type: str = "lgx"
    sort: bool = True

    # --- py-ard DB pinning -------------------------------------------------
    imgt_version: str | None = None  # None = latest available locally

    # --- Error handling ----------------------------------------------------
    #: ``strict`` aborts the batch on the first invalid record (nonzero exit);
    #: ``permissive`` isolates per-record errors and continues.
    mode: str = "permissive"

    # --- Output ------------------------------------------------------------
    output_delimiter: str = "\t"

    def __post_init__(self) -> None:
        if self.redux_type not in VALID_REDUX_TYPES:
            raise ValueError(
                f"invalid redux_type {self.redux_type!r}; "
                f"valid values are {VALID_REDUX_TYPES}"
            )
        if self.mode not in ("strict", "permissive"):
            raise ValueError(f"invalid mode {self.mode!r}; use 'strict' or 'permissive'")
        if self.input_format not in ("auto", "tsv", "csv", "lines"):
            raise ValueError(f"invalid input_format {self.input_format!r}")
        if not self.loci:
            raise ValueError("loci must be non-empty")

    @property
    def qualified_loci(self) -> tuple[str, ...]:
        """Loci with the configured prefix applied (e.g. ``HLA-A``)."""
        return tuple(f"{self.locus_prefix}{locus}" for locus in self.loci)

    def as_manifest_dict(self) -> dict:
        """Serializable view for the provenance manifest."""
        return {
            "gl_column": self.gl_column,
            "donor_column": self.donor_column,
            "input_format": self.input_format,
            "loci": list(self.loci),
            "locus_prefix": self.locus_prefix,
            "enforce_loci": self.enforce_loci,
            "normalize_nomenclature": self.normalize_nomenclature,
            "reduce_alleles": self.reduce_alleles,
            "redux_type": self.redux_type,
            "sort": self.sort,
            "imgt_version": self.imgt_version,
            "mode": self.mode,
        }


@dataclass
class RecordResult:
    """Outcome of normalizing a single input record."""

    donor_id: str | None
    original_gl: str
    normalized_gl: str
    status: str  # "ok" | "error"
    error_detail: str = ""
    loci: tuple[str, ...] = field(default_factory=tuple)
