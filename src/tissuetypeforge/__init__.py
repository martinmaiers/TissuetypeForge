"""TissueTypeForge: simulate HLA typing observations from known genetic truth.

Milestone 1 — observation framework.

Public API:
    parse_gl_string / serialize       GL String round-trip
    Strategy / load_strategy          typing strategy configuration
    observe                           run an observation over truth
    ObservationResult                 structured output
"""

from tissuetypeforge.glstring import parse_gl_string, serialize
from tissuetypeforge.observation import observe
from tissuetypeforge.result import ErrorCategory, ObservationResult
from tissuetypeforge.strategy import Strategy, load_strategy, load_strategy_file

__all__ = [
    "ErrorCategory",
    "ObservationResult",
    "Strategy",
    "load_strategy",
    "load_strategy_file",
    "observe",
    "parse_gl_string",
    "serialize",
]

__version__ = "0.1.0"
