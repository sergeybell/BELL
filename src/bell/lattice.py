from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence

import numpy as np

from .particle import Particle
from .elements import Element


@dataclass(slots=True)
class Lattice:
    elements: List[Element]
    name: str = "LATTICE"

    @classmethod
    def from_elements(cls, elements: Iterable[Element], name: str = "LATTICE") -> "Lattice":
        return cls(list(elements), name=name)

    def total_length(self) -> float:
        return float(sum(getattr(e, "length", 0.0) for e in self.elements))

    def one_turn_map6(self) -> np.ndarray:
        M = np.eye(6, dtype=float)
        for e in self.elements:
            M = e.map6() @ M
        return M

    def track(self, particle: Particle, *, nturns: int = 1, keep_history: bool = True) -> np.ndarray:
        """
        Returns history array of shape (N, 6) if keep_history else final (6,).
        History includes the initial state as the first row.
        """
        x = particle.state.astype(float, copy=True)
        if not keep_history:
            for _ in range(nturns):
                for e in self.elements:
                    x = e.map6() @ x
            return x

        hist: List[np.ndarray] = [x.copy()]
        for _ in range(nturns):
            for e in self.elements:
                x = e.map6() @ x
                hist.append(x.copy())
        return np.vstack(hist)

    def track_through(self, particle: Particle) -> Sequence[np.ndarray]:
        """
        Convenience: track through elements once, returning per-element states (including initial).
        """
        x = particle.state.astype(float, copy=True)
        out = [x.copy()]
        for e in self.elements:
            x = e.map6() @ x
            out.append(x.copy())
        return out

