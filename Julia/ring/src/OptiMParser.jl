# =============================================================================
# OptiMParser.jl - read an OptiM .opt lattice file into a Lattice
# =============================================================================
#
# WHAT
#   parse_optim(path; particle) reads
#     * the "$NAME = expr; => value" variable lines (the value after "=>" is
#       the number OptiM itself evaluated; lines without "=>" must be plain
#       numbers, e.g. "$QUAD =0.1935"),
#     * the header "Energy[MeV]=...  Mass[MeV]=..." and the "Initial:" block
#       (OptiM's initial Twiss - used ONLY as a reference to compare with),
#     * "begin lattice. Number of periods=N" ... "end lattice" (element names),
#     * "begin list" ... "end list" (element parameters "Key[unit]=value").
#   Element type follows the OptiM naming rule (first letter):
#     O... drift, B... bend (field B[kG], gradient G[kG/cm]), Q... quadrupole
#     (G[kG/cm]), S... sextupole (S[kG/cm^2]).  Units are converted to SI and
#     to normalised strengths with B rho from the header energy and mass:
#        angle = L B/(B rho),  k1 = G/(B rho),  k2 = S/(B rho)
#     (S is taken as d2By/dx2 - irrelevant for magnetic_2 where all S = 0).
#
# WHY
#   The program must read the same lattice the OptiM/COSY studies use, with
#   no hand transcription.
#
# NOTE ON "Number of periods"
#   If the listed sequence already bends by 2*pi, it is used ONCE and a
#   warning is stored (magnetic_2.opt says periods=2 but its listed
#   sequence already contains both arcs = 360 deg).  If the listed sequence
#   bends by 2*pi/N it is repeated N times.
#
# WHAT YOU SEE
#   Nothing printed by the parser itself; result.warnings is printed by the
#   scripts.
# =============================================================================
module OptiMParser

using ..Constants, ..Elements

export OptiMFile, parse_optim

struct OptiMFile
    path::String
    variables::Dict{String,Float64}
    energy_MeV::Float64
    mass_MeV::Float64
    initial::Dict{String,Float64}
    nperiods_file::Int
    nrepeat_used::Int
    sequence::Vector{String}
    definitions::Dict{String,Dict{String,Float64}}
    lattice::Lattice
    warnings::Vector{String}
end

const NUM = r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$"

function _value(s::AbstractString, vars, warnings)
    s = strip(s)
    occursin(NUM, s) && return parse(Float64, s)
    sign = 1.0
    if startswith(s, "-")
        sign = -1.0; s = s[2:end]
    end
    if startswith(s, "\$") && haskey(vars, s[2:end])
        return sign * vars[s[2:end]]
    end
    push!(warnings, "cannot evaluate value '$s' - set to 0")
    return 0.0
end

strip_comment(l) = (i = findfirst('#', l); i === nothing ? l : l[1:prevind(l, i)])

function parse_optim(path::AbstractString; particle::Particle = DEUTERON, name = splitext(basename(path))[1])
    lines = readlines(path)
    warnings = String[]
    vars = Dict{String,Float64}()
    i_optim = findfirst(l -> strip(l) == "OptiM", lines)
    i_optim === nothing && error("no 'OptiM' line found in $path")

    # ---- variables (only before the OptiM line) -----------------------------
    for l in lines[1:i_optim-1]
        m = match(r"^\s*\$(\w+)\s*=([^;]*);?\s*(=>\s*(\S+))?", l)
        m === nothing && continue
        if m.captures[4] !== nothing && occursin(NUM, m.captures[4])
            vars[m.captures[1]] = parse(Float64, m.captures[4])
        else
            v = strip(m.captures[2])
            if occursin(NUM, v)
                vars[m.captures[1]] = parse(Float64, v)
            else
                push!(warnings, "variable \$$(m.captures[1]) has no evaluated value")
            end
        end
    end

    # ---- header -------------------------------------------------------------
    text_after = join(lines[i_optim:end], "\n")
    energy = parse(Float64, match(r"Energy\[MeV\]\s*=\s*([-\d.eE+]+)", text_after).captures[1])
    mass   = parse(Float64, match(r"Mass\[MeV\]\s*=\s*([-\d.eE+]+)", text_after).captures[1])
    i_beg = findfirst(l -> occursin("begin lattice", l), lines)
    i_end = findfirst(l -> occursin("end lattice", l), lines)
    initial = Dict{String,Float64}()
    for l in lines[i_optim+1:i_beg-1], m in eachmatch(r"([A-Za-z_]+(\[[^\]]*\])?)\s*=\s*([-+\d.eE]+)", l)
        initial[m.captures[1]] = parse(Float64, m.captures[3])
    end
    mp = match(r"Number of periods\s*=\s*(\d+)", lines[i_beg])
    nper = mp === nothing ? 1 : parse(Int, mp.captures[1])

    # ---- sequence -----------------------------------------------------------
    seq = String[]
    for l in lines[i_beg+1:i_end-1]
        append!(seq, split(strip_comment(l)))
    end

    # ---- element definitions ------------------------------------------------
    i_lb = findfirst(l -> startswith(strip(l), "begin list"), lines)
    i_le = findfirst(l -> startswith(strip(l), "end list"), lines)
    defs = Dict{String,Dict{String,Float64}}()
    for l in lines[i_lb+1:i_le-1]
        t = strip(strip_comment(l))
        isempty(t) && continue
        nm = split(t)[1]
        d = Dict{String,Float64}()
        for m in eachmatch(r"(\w+)\[([^\]]*)\]\s*=\s*(\S+)", t)
            d[m.captures[1]] = _value(m.captures[3], vars, warnings)
        end
        defs[nm] = d
    end

    # Design choice (v0.5): use the physical particle mass (CODATA deuteron), not the
    # OptiM header Mass (2*938.272 MeV rounded).  of.mass_MeV still records the file value.
    ref = reference_from_kinetic(particle, energy)
    brho = ref.brho_Tm
    function make(nm)
        haskey(defs, nm) || error("element $nm used in the sequence but not defined")
        d = defs[nm]
        L = get(d, "L", 0.0) / 100                       # cm -> m
        get(d, "Tilt", 0.0) != 0 && push!(warnings, "$nm: Tilt ignored (not supported in phase 1)")
        c = uppercase(nm[1])
        if c == 'O'
            return Drift(nm, L)
        elseif c == 'B'
            B = get(d, "B", 0.0) / 10                    # kG -> T
            G = get(d, "G", 0.0) * 10                    # kG/cm -> T/m
            B == 0 && G == 0 && return Drift(nm, L)
            return SBend(nm, L, L * B / brho, G / brho)
        elseif c == 'Q'
            return Quadrupole(nm, L, get(d, "G", 0.0) * 10 / brho)
        elseif c == 'S'
            return Sextupole(nm, L, get(d, "S", 0.0) * 1000 / brho)   # kG/cm^2 -> T/m^2
        else
            push!(warnings, "$nm: unknown element type, treated as drift")
            return Drift(nm, L)
        end
    end
    cache = Dict(nm => make(nm) for nm in unique(seq))
    once = [cache[nm] for nm in seq]
    ang = sum(e -> e isa SBend ? e.angle : 0.0, once)
    nrep = 1
    if abs(ang - 2pi) < 1e-6
        nper != 1 && push!(warnings, "file says 'Number of periods=$nper' but the listed sequence already bends by " *
                           "$(round(rad2deg(ang), digits=6)) deg -> sequence used ONCE (as in the COSY deck)")
    elseif abs(ang * nper - 2pi) < 1e-6
        nrep = nper
    else
        push!(warnings, "total bend angle $(rad2deg(ang)) deg (x$nper periods) is not 360 deg")
        nrep = nper
    end
    elems = AbstractElement[]
    for _ in 1:nrep
        append!(elems, once)
    end
    lat = Lattice(name, elems, ref, 1)
    return OptiMFile(String(path), vars, energy, mass, initial, nper, nrep, seq, defs, lat, warnings)
end

end # module
