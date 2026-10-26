#!/usr/bin/env python3
"""Shared, ROOT-independent binning for the SIDIS azimuthal workflows."""

from __future__ import annotations

import math
from typing import Sequence


AZIMUTHAL_BINNING_VERSION = "RP_AZIMUTHAL_V2_15PHI"
AZIMUTHAL_PHI_EDGES = tuple(index * 2.0 * math.pi / 15.0 for index in range(16))
AZIMUTHAL_PT_EDGES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)


def canonical_phi(value: float) -> float:
    """Map a finite angle to the half-open analysis interval [0, 2 pi)."""
    return float(value) % (2.0 * math.pi)


def bin_index(value: float, edges: Sequence[float], close_last: bool = False) -> int | None:
    """Return a one-based bin index, optionally including the final upper edge."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    for index, (low, high) in enumerate(zip(edges[:-1], edges[1:]), 1):
        if low <= value < high or (
            close_last and index == len(edges) - 1 and value == high
        ):
            return index
    return None
