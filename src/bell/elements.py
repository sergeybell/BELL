from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Element(Protocol):
    name: str
    length: float

    def map6(self) -> np.ndarray: ...


def _drift_map(L: float) -> np.ndarray:
    M = np.eye(6, dtype=float)
    M[0, 1] = L
    M[2, 3] = L
    # z, delta untouched in this MVP
    return M


def _quad_map(L: float, k1: float) -> np.ndarray:
    """
    Thick-lens quadrupole linear map in uncoupled x/y.

    x-plane: x'' + k1 x = 0
    y-plane: y'' - k1 y = 0
    """

    def foc_block(k: float) -> np.ndarray:
        if abs(k) < 1e-15:
            return np.array([[1.0, L], [0.0, 1.0]], dtype=float)
        if k > 0:
            r = np.sqrt(k)
            c = np.cos(r * L)
            s = np.sin(r * L)
            return np.array([[c, s / r], [-r * s, c]], dtype=float)
        # k < 0 => defocusing uses cosh/sinh
        r = np.sqrt(-k)
        c = np.cosh(r * L)
        s = np.sinh(r * L)
        return np.array([[c, s / r], [r * s, c]], dtype=float)

    M = np.eye(6, dtype=float)
    Mx = foc_block(k1)
    My = foc_block(-k1)
    M[0:2, 0:2] = Mx
    M[2:4, 2:4] = My
    return M


def _bend_map(L: float, h: float) -> np.ndarray:
    """
    Very simplified sector bend map (horizontal plane only).
    h = 1/rho curvature. This is an MVP placeholder (no dispersion yet).
    """

    # Start with drift-like behavior in both planes.
    M = _drift_map(L)
    # Add weak focusing in x for sector bend: x'' + h^2 x = 0 (approx)
    k = h * h
    if k > 0:
        r = np.sqrt(k)
        c = np.cos(r * L)
        s = np.sin(r * L)
        Mx = np.array([[c, s / r], [-r * s, c]], dtype=float)
        M[0:2, 0:2] = Mx
    return M


@dataclass(frozen=True, slots=True)
class Drift:
    length: float
    name: str = "DRIFT"

    def map6(self) -> np.ndarray:
        return _drift_map(self.length)


@dataclass(frozen=True, slots=True)
class Quadrupole:
    length: float
    k1: float
    name: str = "QUAD"

    def map6(self) -> np.ndarray:
        return _quad_map(self.length, self.k1)


@dataclass(frozen=True, slots=True)
class Bend:
    length: float
    h: float
    name: str = "BEND"

    def map6(self) -> np.ndarray:
        return _bend_map(self.length, self.h)

