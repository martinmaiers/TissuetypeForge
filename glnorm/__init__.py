"""glnorm -- a standalone GL String normalization CLI harness.

Independent of TissueTypeForge / DiplotypeForge. Depends only on ``glstring``,
``py-ard``, ``gl_smartsort`` and the standard library.

Three independently toggleable operations:

* nomenclature normalization (lossless, glstring)
* allele reduction (optional, lossy, py-ard)
* canonical ordering (gl_smartsort)
"""

from __future__ import annotations

__version__ = "0.1.0"

from .config import DEFAULT_LOCI, VALID_REDUX_TYPES, NormalizeConfig, RecordResult
from .normalize import (
    detect_loci,
    normalize_nomenclature,
    reduce_alleles,
    sort_canonical,
    validate,
)
from .pipeline import normalize_one, run_pipeline

__all__ = [
    "__version__",
    "NormalizeConfig",
    "RecordResult",
    "DEFAULT_LOCI",
    "VALID_REDUX_TYPES",
    "validate",
    "detect_loci",
    "normalize_nomenclature",
    "reduce_alleles",
    "sort_canonical",
    "normalize_one",
    "run_pipeline",
]
