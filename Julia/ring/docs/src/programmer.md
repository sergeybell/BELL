# Programmer manual (stub)

To be written. Entry point: `src/LieRing.jl` includes all modules in
dependency order and re-exports their names; scripts do
`include("../src/LieRing.jl"); using .LieRing`. Build these docs with
`julia --project=docs docs/make.jl` (HTML in `docs/build/`) and the PDF with
`julia docs/make_pdf.jl` (needs xelatex).
