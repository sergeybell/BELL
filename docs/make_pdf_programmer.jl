# Programmer Manual PDF
empty!(ARGS)
push!(ARGS, "programmer")
include(joinpath(@__DIR__, "make_pdf.jl"))
