# Spin

Spin is not part of the present equations of motion.

The species record stores the anomalous magnetic moment ``G`` next to the mass and the charge. A particle state carries a spin slot beside the canonical 6-vector. The slot is empty unless a caller fills it, and every tracker leaves it untouched.

The equation that belongs here is the Thomas–Bargmann–Michel–Telegdi equation in the same Frenet frame as ``z``, with the laboratory field taken from the element that already defines ``H``. That equation is not evaluated. Orbit, tunes, chromaticity and the RF bucket do not depend on ``G``.

Until the equation exists, ``G`` is only a stored constant, and the spin slot is only storage.
