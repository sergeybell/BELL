# =============================================================================
# optim_converter.jl - write a Bell structure file from an OptiM .opt
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   convert_optim(opt_path, structure_path) reads an OptiM lattice file and
#   writes a Bell structure file.  It does not assemble a Lattice: that is
#   load_structure's job.
#   Element type follows the OptiM naming rule (first letter):
#     O... drift, B... bend, Q... quadrupole, S... sextupole.
#   Lengths cm -> m, fields kG -> T, gradients to normalised k1, k2 using
#   B rho from the header kinetic energy and the CODATA particle mass.
#
# =============================================================================
module OptimConverter

using ..Coordinates, ..Elements, ..Structure

export OptiMFile, convert_optim

"""
    OptiMFile

Result of converting an OptiM `.opt` file ([`convert_optim`](@ref)).

# Fields
- `path::String` – the source `.opt` file; `structure_path::String` – written `.bell` file.
- `variables::Dict{String,Float64}` – OptiM variables (dollar-prefixed names) after evaluation.
- `energy_MeV::Float64`, `mass_MeV::Float64` – kinetic energy and mass from the header.
- `initial::Dict{String,Float64}` – initial Twiss values given in the file.
- `nperiods_file::Int` – number of periods declared; `nrepeat_used::Int` – repetitions written to the ring.
- `sequence::Vector{String}` – element names of one period.
- `definitions::Dict{String,Dict{String,Float64}}` – raw element parameters (OptiM units).
- `warnings::Vector{String}` – everything that was ignored or guessed.

# Conversion rules
The first letter of an element name selects the kind: O = drift, B = sector bend,
Q = quadrupole, S = sextupole (other letters are reported in `warnings`).
Units: L [cm] → m (/100); bend field B [kG] → T (/10), angle = L·B/(Bρ);
gradient G [kG/cm] → T/m (×10), k1 = G/(Bρ); sextupole S [kG/cm²] → T/m² (×1000),
k2 = S/(Bρ). A B-element with B = G = 0 becomes a drift; `Tilt` is ignored (warning).
Bρ is computed from the header energy and the CODATA mass of `particle`
(the header `Mass` is recorded but not used). Periods: if the listed sequence
already bends by 360° it is used once, else it is repeated `Number of periods`
times when that closes the ring (otherwise used once with a warning).
"""
struct OptiMFile
    path::String
    structure_path::String
    variables::Dict{String,Float64}
    energy_MeV::Float64
    mass_MeV::Float64
    initial::Dict{String,Float64}
    nperiods_file::Int
    nrepeat_used::Int
    sequence::Vector{String}
    definitions::Dict{String,Dict{String,Float64}}
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

function _parse_optim_raw(path::AbstractString; particle::Particle = DEUTERON)
    lines = readlines(path)
    warnings = String[]
    vars = Dict{String,Float64}()
    i_optim = findfirst(l -> strip(l) == "OptiM", lines)
    i_optim === nothing && error("no 'OptiM' line found in $path")

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

    seq = String[]
    for l in lines[i_beg+1:i_end-1]
        append!(seq, split(strip_comment(l)))
    end

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

    ref = reference_from_kinetic(particle, energy)
    brho = ref.brho_Tm
    function make(nm)
        haskey(defs, nm) || error("element $nm used in the sequence but not defined")
        d = defs[nm]
        L = get(d, "L", 0.0) / 100
        get(d, "Tilt", 0.0) != 0 && push!(warnings, "$nm: Tilt ignored")
        c = uppercase(nm[1])
        if c == 'O'
            return Drift(nm, L)
        elseif c == 'B'
            B = get(d, "B", 0.0) / 10
            G = get(d, "G", 0.0) * 10
            B == 0 && G == 0 && return Drift(nm, L)
            return SBend(nm, L, L * B / brho, G / brho)
        elseif c == 'Q'
            return Quadrupole(nm, L, get(d, "G", 0.0) * 10 / brho)
        elseif c == 'S'
            return Sextupole(nm, L, get(d, "S", 0.0) * 1000 / brho)
        else
            push!(warnings, "$nm: unknown element type, treated as drift")
            return Drift(nm, L)
        end
    end
    cache = Dict(nm => make(nm) for nm in unique(seq))
    once = [cache[nm] for nm in seq]
    ang = sum(e -> e isa SBend ? e.angle : 0.0, once)
    nrep = 1
    if abs(ang - 2pi) < deg2rad(1.0)
        nper != 1 && push!(warnings, "file says 'Number of periods=$nper' but the listed sequence already bends by " *
                           "$(round(rad2deg(ang), digits=6)) deg -> sequence used ONCE")
    elseif abs(ang * nper - 2pi) < deg2rad(1.0)
        nrep = nper
        push!(warnings, "sequence bend $(round(rad2deg(ang),digits=4)) deg repeated nper=$nper to make 360 deg")
    else
        push!(warnings, "total bend angle $(rad2deg(ang)) deg (x$nper periods) is not 360 deg; using sequence once")
        nrep = 1
    end
    elems = AbstractElement[]
    for _ in 1:nrep
        append!(elems, once)
    end
    lat = Lattice(splitext(basename(String(path)))[1], elems, ref, 1)
    meta = (path = String(path), variables = vars, energy_MeV = energy, mass_MeV = mass,
            initial = initial, nperiods_file = nper, nrepeat_used = nrep, sequence = seq,
            definitions = defs, warnings = warnings, lattice = lat)
    return meta
end

"""
    convert_optim(opt_path, structure_path; particle = DEUTERON) -> OptiMFile

Read an OptiM lattice file, write an equivalent Bell structure file to
`structure_path`, and return conversion metadata ([`OptiMFile`](@ref)).
The `Lattice` is obtained afterwards with `load_structure(structure_path)`;
Bell itself never reads OptiM files at run time.
"""
function convert_optim(opt_path::AbstractString, structure_path::AbstractString;
                       particle::Particle = DEUTERON)
    meta = _parse_optim_raw(opt_path; particle)
    extra = "# converted from $(meta.path)\n# OptiM header Mass = $(meta.mass_MeV) MeV (not used; CODATA particle mass)"
    write_structure(structure_path, meta.lattice; extra_header = extra)
    return OptiMFile(meta.path, String(structure_path), meta.variables, meta.energy_MeV, meta.mass_MeV,
                     meta.initial, meta.nperiods_file, meta.nrepeat_used, meta.sequence,
                     meta.definitions, meta.warnings)
end

end # module
