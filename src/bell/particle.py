from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class Particle:
    """
    6D canonical-ish coordinates vector:
      (x, px, y, py, z, delta)

    This MVP uses a linear map approach; exact conventions can evolve later.
    """

    state: np.ndarray  # shape (6,)

    @classmethod
    def from_coords(
        cls,
        x: float = 0.0,
        px: float = 0.0,
        y: float = 0.0,
        py: float = 0.0,
        z: float = 0.0,
        delta: float = 0.0,
    ) -> "Particle":
        return cls(np.array([x, px, y, py, z, delta], dtype=float))

    def copy(self) -> "Particle":
        return Particle(self.state.copy())

