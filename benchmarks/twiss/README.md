# Twiss cross-code benchmark — magnetic_2

Same lattice (examples/lattices/magnetic_2.opt, deuteron, E_kin = 270 MeV, C = 117.388 m) in four codes.
Each subfolder: `run.jl` (runs the code headless, parses output -> `results.csv`, `timings.csv`), `raw/`, `README.md`.

    julia --project=examples benchmarks/twiss/optim/run.jl   # each run.jl also writes raw/twiss_std.csv + twiss_<code>.png
    julia --project=examples benchmarks/twiss/bell/run.jl
    julia --project=examples benchmarks/twiss/cosy/run.jl
    julia --project=examples benchmarks/twiss/madx/run.jl
    julia --project=examples benchmarks/twiss/compare.jl            # -> comparison.md
    julia --project=examples benchmarks/twiss/plot_compare.jl       # -> twiss_overlay.png, twiss_diff_vs_bell*.png, twiss_diff_summary.csv
    julia --project=examples benchmarks/twiss/bell/linear_routes.jl all   # Bell: ring_matrix vs Dragt-Finn vs LieFlow N=1/2

`common.jl` — shared CSV/plot helpers (Plots.jl from examples env; loading Plots adds ~3-4 s to each script, not part of code timings).

**`twiss_benchmark.ipynb` is TEMPORARY** (review notebook, Julia kernel, reads the files above; to be removed later).

## (a) What each code computes here (headless)

| quantity | OptiM CLI | Bell | COSY 10.2 | MAD-X 5 |
|---|---|---|---|---|
| Twiss beta/alpha along ring | yes (-f/-s) | yes | yes (beta, D) | yes |
| Dispersion D, D' | yes, but from Initial (not periodic) | yes, periodic | yes | yes |
| Tunes Q | yes (accumulated phase) | yes | fractional (TP) | yes |
| Chromaticity xi | GUI only (View/Integrals) | yes (order-3 map) | yes (TP, OV 3 2 1) | yes (TWISS CHROM, PTC) |
| alpha_c / eta | GUI only | yes | yes (from 6D map) | yes |
| Radiation integrals I1..I5 | GUI only | from Twiss (script) | not scripted | yes (SYNCH_1..5) |
| Higher order (D2, xi2, alpha1) | no | yes | yes (higher OV) | PTC |

## (b) Values: see comparison.md (generated).

## (c) Timings: see comparison.md; pure computation vs whole process is summarised in the top README of the report.

OptiM GUI numbers given by the user are kept only in optim/gui_reference_manual.txt (manual reference, not a run).
