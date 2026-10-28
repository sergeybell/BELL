# =============================================================================
# structure.jl - structure file, Lattice, and lattice-editing helpers
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   A structure file has three parts:
#     reference    particle, kinetic energy, harmonic
#     parameters   named numbers (the optimizer changes these)
#     elements     one card per named magnet / drift / cavity
#     sequence     names only; `N * cell` repeats a named sequence
#   load_structure(path) reads the file into a Lattice.
#   write_structure(path, lat) dumps a Lattice back to the same format.
#   Lattice-editing helpers (sextupole families, cavity insert) live here:
#   they change the sequence, not the element types.
#
# =============================================================================
module Structure

using Printf
using ..Coordinates, ..Elements

export Lattice, circumference, total_bend_angle, unique_elements,
       SEXT_FAMILIES, find_optics_maxima, place_sextupole_families,
       set_family_strengths, insert_cavity, set_cavity_phase, with_reference,
       load_structure, write_structure, StructureFile

"A ring: ordered elements, reference particle, harmonic number used for phi."
struct Lattice
    name::String
    elements::Vector{AbstractElement}
    ref::RefParticle
    harmonic::Int
end

circumference(lat::Lattice) = sum(elem_length, lat.elements)
total_bend_angle(lat::Lattice) = sum(e -> e isa SBend ? e.angle : 0.0, lat.elements)

"Distinct elements (by value) in order of first appearance."
function unique_elements(lat::Lattice)
    seen = AbstractElement[]
    for e in lat.elements
        any(isequal(e), seen) || push!(seen, e)
    end
    return seen
end

const SEXT_FAMILIES = ["SXD", "SXB", "SYB"]

"""
find_optics_maxima(opt; frac = 0.8, dmin = 0.5)
Local maxima of D_x, beta_x, beta_y in the optics table `opt.data`.
"""
function find_optics_maxima(opt; frac = 0.8, dmin = 0.5)
    d = opt.data
    found = NamedTuple[]
    for (q, arr) in (("D_x", d.D), ("beta_x", d.betx), ("beta_y", d.bety))
        mask = d.D .> dmin
        g = maximum(arr[mask])
        last_s = -1e9
        for k in 2:length(arr)-1
            (mask[k] && arr[k] >= arr[k-1] && arr[k] >= arr[k+1] && arr[k] >= frac * g) || continue
            d.s[k] - last_s < 0.5 && continue
            last_s = d.s[k]
            push!(found, (quantity = q, s = d.s[k], value = arr[k], betx = d.betx[k], bety = d.bety[k], D = d.D[k]))
        end
    end
    return found
end

"""
place_sextupole_families(lat, opt) -> (new_lattice, fam_idx, rows, notes)
Puts a sextupole of family SXD / SXB / SYB (strength 0) into the nearest
short slot of each D_x / beta_x / beta_y maximum.
"""
function place_sextupole_families(lat::Lattice, opt; kw...)
    maxima = find_optics_maxima(opt; kw...)
    fam_of = Dict("D_x" => "SXD", "beta_x" => "SXB", "beta_y" => "SYB")
    assigned = Dict{Int,String}()
    notes = String[]; rows = NamedTuple[]
    for q in ("D_x", "beta_x", "beta_y"), m in filter(m -> m.quantity == q, maxima)
        best, bd = 0, 1e9
        for (i, (s0, s1, e)) in enumerate(opt.bounds)
            (e isa Sextupole || (e isa Drift && e.L <= 0.5)) || continue
            dist = abs((s0 + s1) / 2 - m.s)
            dist < bd && ((best, bd) = (i, dist))
        end
        if best == 0 || bd > 1.0
            push!(notes, "no slot within 1 m of $(q) max at s=$(round(m.s, digits=3))"); continue
        end
        if haskey(assigned, best)
            push!(notes, "$(q) max at s=$(round(m.s, digits=3)) m shares slot $(opt.bounds[best][3].name) with family $(assigned[best])")
            continue
        end
        assigned[best] = fam_of[q]
        push!(rows, (family = fam_of[q], quantity = q, s_max = m.s, value = m.value, slot = opt.bounds[best][3].name,
                     s_slot = (opt.bounds[best][1] + opt.bounds[best][2]) / 2, betx = m.betx, bety = m.bety, D = m.D))
    end
    elems = AbstractElement[]
    fam_idx = Dict(f => Int[] for f in SEXT_FAMILIES)
    for (i, e) in enumerate(lat.elements)
        if haskey(assigned, i)
            f = assigned[i]; L = elem_length(e)
            if L <= 0.15
                push!(elems, Sextupole(f, L, 0.0)); push!(fam_idx[f], length(elems))
            else
                push!(elems, Drift(e.name * "a", (L - 0.1) / 2)); push!(elems, Sextupole(f, 0.1, 0.0))
                push!(fam_idx[f], length(elems)); push!(elems, Drift(e.name * "b", (L - 0.1) / 2))
            end
        else
            push!(elems, e)
        end
    end
    return Lattice(lat.name * "+sext", elems, lat.ref, lat.harmonic), fam_idx, rows, notes
end

"New lattice with sextupole strengths k2[family] (families not given -> 0)."
function set_family_strengths(lat::Lattice, fam_idx, k2::Dict)
    elems = copy(lat.elements)
    for (f, idx) in fam_idx, i in idx
        elems[i] = Sextupole(f, elems[i].L, get(k2, f, 0.0))
    end
    return Lattice(lat.name, elems, lat.ref, lat.harmonic)
end

"Split drift `in_element` (first occurrence) as drift - cavity - drift."
function insert_cavity(lat::Lattice, cav::Cavity; in_element = "OBM")
    i = findfirst(e -> e.name == in_element, lat.elements)
    e = lat.elements[i]
    Ld = (e.L - cav.L) / 2
    elems = vcat(lat.elements[1:i-1], [Drift(e.name * "_u", Ld), cav, Drift(e.name * "_d", Ld)], lat.elements[i+1:end])
    return Lattice(lat.name * "+RF", elems, lat.ref, cav.harmonic), i + 1
end

"Same lattice with every finite cavity set to synchronous phase phi_s."
set_cavity_phase(lat::Lattice, phi_s) =
    Lattice(lat.name, [e isa Cavity ? Cavity(e.name, e.L, e.voltage_MV, e.harmonic, phi_s) : e for e in lat.elements],
            lat.ref, lat.harmonic)

"Same elements with a different reference particle."
with_reference(lat::Lattice, ref::RefParticle) = Lattice(lat.name, lat.elements, ref, lat.harmonic)

# =============================================================================
# Structure-file parser
# =============================================================================

struct StructureFile
    path::String
    parameters::Dict{String,Float64}
    lattice::Lattice
end

_strip_comment(l) = (i = findfirst('#', l); i === nothing ? String(l) : String(l[1:prevind(l, i)]))

const _PARTICLES = Dict(
    "deuteron" => DEUTERON,
    "proton"   => Particle("proton", 938.27208816, 1.0, 1.7928473443),
    "electron" => Particle("electron", 0.51099895000, -1.0, 0.001159652181),
)

function _tok(s::AbstractString)
    s = strip(s)
    isempty(s) && return String[]
    return String.(split(s))
end

"Evaluate a numeric expression against already-defined parameters."
function _eval_expr(expr::AbstractString, params::Dict{String,Float64})
    s = strip(expr)
    isempty(s) && error("empty expression")
    deg = false
    if endswith(lowercase(s), "deg")
        deg = true
        s = strip(s[1:end-3])
    end
    val = _eval_sum(s, params)
    return deg ? deg2rad(val) : val
end

function _eval_sum(s::AbstractString, params)
    parts, ops = String[], Char[]
    depth = 0
    buf = IOBuffer()
    for c in s
        if c == '('
            depth += 1; write(buf, c)
        elseif c == ')'
            depth -= 1; write(buf, c)
        elseif depth == 0 && (c == '+' || c == '-') && position(buf) > 0
            push!(parts, String(take!(buf)))
            push!(ops, c)
        else
            write(buf, c)
        end
    end
    push!(parts, String(take!(buf)))
    val = _eval_prod(strip(parts[1]), params)
    for (op, p) in zip(ops, parts[2:end])
        v = _eval_prod(strip(p), params)
        val = op == '+' ? val + v : val - v
    end
    return val
end

function _eval_prod(s::AbstractString, params)
    s = strip(s)
    parts, ops = String[], Char[]
    depth = 0
    buf = IOBuffer()
    for c in s
        if c == '('
            depth += 1; write(buf, c)
        elseif c == ')'
            depth -= 1; write(buf, c)
        elseif depth == 0 && (c == '*' || c == '/')
            push!(parts, String(take!(buf)))
            push!(ops, c)
        else
            write(buf, c)
        end
    end
    push!(parts, String(take!(buf)))
    val = _eval_atom(strip(parts[1]), params)
    for (op, p) in zip(ops, parts[2:end])
        v = _eval_atom(strip(p), params)
        val = op == '*' ? val * v : val / v
    end
    return val
end

function _eval_atom(s::AbstractString, params)
    s = strip(s)
    isempty(s) && error("empty atom in expression")
    if startswith(s, "(") && endswith(s, ")")
        return _eval_sum(s[2:end-1], params)
    end
    if startswith(s, "-") && length(s) > 1
        return -_eval_atom(s[2:end], params)
    end
    if startswith(s, "+") && length(s) > 1
        return _eval_atom(s[2:end], params)
    end
    if occursin(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$", s)
        return parse(Float64, s)
    end
    haskey(params, s) || error("unknown parameter '$s'")
    return params[s]
end

function _parse_kv(tokens, i, params)
    kv = Dict{String,Float64}()
    while i <= length(tokens)
        t = tokens[i]
        if i < length(tokens) && tokens[i+1] == "="
            j = i + 2
            buf = String[]
            while j <= length(tokens)
                if j + 1 <= length(tokens) && tokens[j+1] == "="
                    break
                end
                push!(buf, tokens[j]); j += 1
            end
            isempty(buf) && error("missing value for $(tokens[i])")
            kv[lowercase(t)] = _eval_expr(join(buf, " "), params)
            i = j
        else
            error("expected key = value, got '$(tokens[i])'")
        end
    end
    return kv
end

function _make_element(name, kind, kv)
    k = lowercase(kind)
    L = get(kv, "l", 0.0)
    if k == "drift"
        return Drift(name, L)
    elseif k == "sbend" || k == "bend"
        return SBend(name, L, get(kv, "angle", 0.0), get(kv, "k1", 0.0))
    elseif k == "quadrupole" || k == "quad"
        return Quadrupole(name, L, get(kv, "k1", 0.0))
    elseif k == "sextupole" || k == "sext"
        return Sextupole(name, L, get(kv, "k2", 0.0))
    elseif k == "rfcavity" || k == "rf"
        V = haskey(kv, "v") ? kv["v"] : get(kv, "voltage_mv", 0.0)
        h = Int(round(haskey(kv, "h") ? kv["h"] : get(kv, "harmonic", 1.0)))
        ph = get(kv, "phi_s", get(kv, "phis", 0.0))
        return RFCavity(name, V, h, ph)
    elseif k == "cavity"
        V = haskey(kv, "v") ? kv["v"] : get(kv, "voltage_mv", 0.0)
        h = Int(round(haskey(kv, "h") ? kv["h"] : get(kv, "harmonic", 1.0)))
        ph = get(kv, "phi_s", get(kv, "phis", 0.0))
        return Cavity(name, L, V, h, ph)
    elseif k == "marker"
        return Marker(name)
    else
        error("unknown element kind '$kind' for $name")
    end
end

function _expand_seq(tokens, elements, sequences; stack = String[])
    out = AbstractElement[]
    i = 1
    while i <= length(tokens)
        t = tokens[i]
        if i + 1 <= length(tokens) && tokens[i+1] == "*"
            n = Int(round(parse(Float64, t)))
            i + 2 <= length(tokens) || error("expected name after '*'")
            name = tokens[i+2]
            piece = _expand_name(name, elements, sequences, stack)
            for _ in 1:n
                append!(out, piece)
            end
            i += 3
        elseif t == "*"
            error("'*' without a count")
        else
            append!(out, _expand_name(t, elements, sequences, stack))
            i += 1
        end
    end
    return out
end

function _expand_name(name, elements, sequences, stack)
    name in stack && error("cyclic sequence '$name'")
    if haskey(elements, name)
        return AbstractElement[elements[name]]
    elseif haskey(sequences, name)
        return _expand_seq(sequences[name], elements, sequences; stack = vcat(stack, name))
    else
        error("name '$name' is not an element and not a sequence")
    end
end

"""
    load_structure(path) -> Lattice

Read a Bell structure file.
"""
function load_structure(path::AbstractString)
    lines = readlines(path)
    section = ""
    seq_name = ""
    ref_kv = Dict{String,String}()
    params = Dict{String,Float64}()
    elements = Dict{String,AbstractElement}()
    elem_order = String[]
    sequences = Dict{String,Vector{String}}()
    ring_tokens = String[]
    pending_seq = String[]

    function flush_seq!()
        isempty(pending_seq) && return
        if isempty(seq_name)
            append!(ring_tokens, pending_seq)
        else
            sequences[seq_name] = copy(pending_seq)
        end
        empty!(pending_seq)
    end

    for raw in lines
        line = strip(_strip_comment(raw))
        isempty(line) && continue
        low = lowercase(line)
        if low == "reference" || startswith(low, "reference ")
            flush_seq!(); section = "reference"; continue
        elseif low == "parameters" || startswith(low, "parameters ")
            flush_seq!(); section = "parameters"; continue
        elseif low == "elements" || startswith(low, "elements ")
            flush_seq!(); section = "elements"; continue
        elseif startswith(low, "sequence")
            flush_seq!()
            rest = strip(line[length("sequence")+1:end])
            seq_name = rest
            section = "sequence"
            continue
        end
        if section == "reference"
            toks = _tok(line)
            length(toks) >= 2 || error("reference line needs a key and a value: $line")
            # `particle deuteron` or `kinetic_MeV 270`
            ref_kv[lowercase(toks[1])] = join(toks[2:end], " ")
        elseif section == "parameters"
            m = match(r"^([A-Za-z_][\w]*)\s*=\s*(.+)$", line)
            m === nothing && error("parameter line must be NAME = EXPR: $line")
            params[m.captures[1]] = _eval_expr(m.captures[2], params)
        elseif section == "elements"
            toks = _tok(line)
            length(toks) >= 2 || error("element card needs NAME KIND: $line")
            nm, kind = toks[1], toks[2]
            kv = isempty(toks[3:end]) ? Dict{String,Float64}() : _parse_kv(toks, 3, params)
            elements[nm] = _make_element(nm, kind, kv)
            push!(elem_order, nm)
        elseif section == "sequence"
            seqline = replace(line, r"(\d+)\s*\*\s*" => s"\1 * ")
            append!(pending_seq, _tok(seqline))
        else
            error("line outside a section: $line")
        end
    end
    flush_seq!()

    haskey(ref_kv, "particle") || error("structure file missing 'particle'")
    haskey(ref_kv, "kinetic_mev") || error("structure file missing 'kinetic_MeV'")
    pname = lowercase(ref_kv["particle"])
    haskey(_PARTICLES, pname) || error("unknown particle '$pname'")
    particle = _PARTICLES[pname]
    Ekin = _eval_expr(ref_kv["kinetic_mev"], params)
    mass = haskey(ref_kv, "mass_mev") ? _eval_expr(ref_kv["mass_mev"], params) : particle.mass_MeV
    harm = haskey(ref_kv, "harmonic") ? Int(round(_eval_expr(ref_kv["harmonic"], params))) : 1
    ref = reference_from_kinetic(particle, Ekin; mass_MeV = mass)

    isempty(ring_tokens) && error("structure file has no unnamed sequence (the ring)")
    elems = _expand_seq(ring_tokens, elements, sequences)
    name = splitext(basename(String(path)))[1]
    lat = Lattice(name, elems, ref, harm)
    return lat
end

_fmt(x::Real) = @sprintf("%.16g", x)

function _card_line(e::AbstractElement)
    nm = rpad(e.name, 8)
    if e isa Drift
        return @sprintf("%s %-12s L = %-20s", nm, "drift", _fmt(e.L))
    elseif e isa SBend
        return @sprintf("%s %-12s L = %-20s angle = %-20s k1 = %s",
                        nm, "sbend", _fmt(e.L), _fmt(e.angle), _fmt(e.k1))
    elseif e isa Quadrupole
        return @sprintf("%s %-12s L = %-20s k1 = %s", nm, "quadrupole", _fmt(e.L), _fmt(e.k1))
    elseif e isa Sextupole
        return @sprintf("%s %-12s L = %-20s k2 = %s", nm, "sextupole", _fmt(e.L), _fmt(e.k2))
    elseif e isa RFCavity
        return @sprintf("%s %-12s V = %-20s h = %-6d phi_s = %s",
                        nm, "rfcavity", _fmt(e.voltage_MV), e.harmonic, _fmt(e.phi_s))
    elseif e isa Cavity
        return @sprintf("%s %-12s L = %-20s V = %-20s h = %-6d phi_s = %s",
                        nm, "cavity", _fmt(e.L), _fmt(e.voltage_MV), e.harmonic, _fmt(e.phi_s))
    elseif e isa Marker
        return @sprintf("%-8s marker", e.name)
    else
        error("cannot write element $(typeof(e))")
    end
end

"""
    write_structure(path, lat; parameters, sequences)

Write a Lattice as a structure file.  Sequence is the flat list of element
names unless `sequences` provides named blocks.
"""
function write_structure(path::AbstractString, lat::Lattice;
                         parameters::Dict{String,Float64} = Dict{String,Float64}(),
                         extra_header = "")
    io = IOBuffer()
    println(io, "# Bell structure file")
    println(io, "# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>")
    isempty(extra_header) || println(io, extra_header)
    println(io)
    println(io, "reference")
    println(io, "  particle ", lat.ref.particle.name)
    println(io, "  kinetic_MeV ", lat.ref.Ekin_MeV)
    println(io, "  harmonic ", lat.harmonic)
    if lat.ref.mass_MeV != lat.ref.particle.mass_MeV
        println(io, "  mass_MeV ", lat.ref.mass_MeV)
    end
    println(io)
    if !isempty(parameters)
        println(io, "parameters")
        for (k, v) in sort(collect(parameters))
            println(io, "  ", k, " = ", v)
        end
        println(io)
    end
    println(io, "elements")
    seen = Dict{String,AbstractElement}()
    names = String[]
    for e in lat.elements
        if !haskey(seen, e.name)
            seen[e.name] = e
            println(io, "  ", _card_line(e))
        elseif !isequal(seen[e.name], e)
            error("element name $(e.name) used with two different parameter sets")
        end
        push!(names, e.name)
    end
    println(io)
    println(io, "sequence")
    # wrap names ~ 8 per line
    for i in 1:8:length(names)
        println(io, "  ", join(names[i:min(i+7, length(names))], " "))
    end
    write(path, String(take!(io)))
    return path
end

end # module
