"""Optional publication-figure QA helpers used by repository plotting scripts.

The full panel-alignment auditor is supplied by the local figure-making skill
when available.  It is intentionally not vendored into this research
repository.  These wrappers keep the plotting scripts portable: a missing
auditor produces a clear warning, while a configured auditor is still run
with the same strict arguments used during manuscript figure generation.
"""

from __future__ import annotations

import importlib
import warnings
from typing import Any


def require_matplotlib_panel_alignment(fig: Any, **kwargs: Any) -> Any:
    """Run the optional panel-alignment gate, if it is installed.

    The repository's numerical output is independent of the external QA
    package.  Callers should still treat a ``None`` return as *not audited*;
    it is never a claim that alignment passed.
    """

    try:
        auditor = importlib.import_module("audit_panel_alignment")
    except ModuleNotFoundError:
        warnings.warn(
            "Optional audit_panel_alignment is not installed; panel alignment "
            "was not audited for this local export.",
            RuntimeWarning,
            stacklevel=2,
        )
        return None
    return auditor.require_matplotlib_panel_alignment(fig, **kwargs)
