from __future__ import annotations

import numpy as np

from bell import Bend, Drift, Lattice, Particle, Quadrupole, twiss_from_2x2


def main() -> None:
    # Simple FODO-like cell with a bend placeholder
    lat = Lattice.from_elements(
        [
            Drift(1.0, name="D1"),
            Quadrupole(0.3, k1=+1.2, name="QF"),
            Drift(0.5, name="D2"),
            Bend(0.8, h=0.2, name="B1"),
            Drift(0.5, name="D3"),
            Quadrupole(0.3, k1=-1.2, name="QD"),
            Drift(1.0, name="D4"),
        ],
        name="demo",
    )

    M = lat.one_turn_map6()
    Mx = M[0:2, 0:2]
    My = M[2:4, 2:4]

    print(f"Lattice length: {lat.total_length():.3f} m")
    print("One-turn 2x2 maps:")
    print("Mx=\n", np.array2string(Mx, precision=6, suppress_small=True))
    print("My=\n", np.array2string(My, precision=6, suppress_small=True))

    try:
        tx = twiss_from_2x2(Mx)
        ty = twiss_from_2x2(My)
        print(f"Twiss X: beta={tx.beta:.6g}, alpha={tx.alpha:.6g}, mu={tx.mu:.6g} rad")
        print(f"Twiss Y: beta={ty.beta:.6g}, alpha={ty.alpha:.6g}, mu={ty.mu:.6g} rad")
    except ValueError as e:
        print("Twiss failed:", e)

    p0 = Particle.from_coords(x=1e-3, px=0.0, y=2e-3, py=0.0)
    hist = lat.track(p0, nturns=3, keep_history=True)
    print("Tracking history (first 6 rows):")
    print(np.array2string(hist[:6, :4], precision=6, suppress_small=True))


if __name__ == "__main__":
    main()
