from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class Twiss:
    beta: float
    alpha: float
    mu: float  # phase advance [rad]

    @property
    def gamma(self) -> float:
        return (1.0 + self.alpha * self.alpha) / self.beta


def twiss_from_2x2(M: np.ndarray) -> Twiss:
    """
    Compute periodic Twiss for a stable 2x2 symplectic-ish matrix M.

    For stable motion: |trace(M)| < 2.
    """
    M = np.asarray(M, dtype=float)
    if M.shape != (2, 2):
        raise ValueError("Expected 2x2 matrix")

    tr = float(M[0, 0] + M[1, 1])
    if abs(tr) >= 2.0:
        raise ValueError(f"Unstable motion: trace={tr:.6g}")

    mu = float(np.arccos(tr / 2.0))
    sin_mu = float(np.sin(mu))
    if abs(sin_mu) < 1e-15:
        raise ValueError("Degenerate phase advance (sin(mu) ~ 0)")

    beta = float(M[0, 1] / sin_mu)
    alpha = float((M[0, 0] - M[1, 1]) / (2.0 * sin_mu))
    if beta <= 0:
        raise ValueError(f"Non-physical beta computed: {beta:.6g}")
    return Twiss(beta=beta, alpha=alpha, mu=mu)

