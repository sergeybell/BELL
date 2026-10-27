# Linear optics

The linear map of an element is ``R=\exp(L\,S\,\mathcal{H})`` with
``\mathcal{H}_{ij}=\partial^2 H/\partial z_i\partial z_j|_0``. The one-turn
matrix ``M`` is the product over the ring. With ``\cos\mu=\tfrac12\mathrm{Tr}M_x``,

```math
\beta=\frac{M_{12}}{\sin\mu},\qquad
\alpha=\frac{M_{11}-M_{22}}{2\sin\mu},
```

and the closed dispersion solves ``(I-M_4)D=M_{\cdot 6}``. Momentum compaction
to first order:

```math
\alpha_0=\frac{1}{\gamma_0^2}-\frac{(l|\delta)_{\mathrm{closed}}}{C}
=\frac1C\oint h\,D\,\mathrm{d}s.
```

Twiss and dispersion propagated in slices close after one turn to machine
precision.

![Closed beta functions, horizontal dispersion and magnet layout.](assets/magnetic_2_twiss.png)