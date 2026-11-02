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
#   write_structure(path, lattice) dumps a Lattice back to the same format.
#   insert_cavity splits a named drift and inserts a thin RFCavity. Sextupole families live in
#   sextupole_families.jl.
#
# =============================================================================
module Structure

using Printf
using ..Coordinates, ..Elements

export Lattice, path_length, total_bend_angle, unique_elements,
       insert_cavity, set_cavity_phase, with_reference,
       load_structure, write_structure, StructureFile

"""
    Lattice(name, elements, ref, harmonic)

A ring: ordered list of elements, reference particle and RF harmonic number.

# Fields
- `name::String` – lattice name (file base name for `load_structure`).
- `elements::Vector{AbstractElement}` – the full ring in tracking order (sequences already expanded; repeated cells are the same element values).
- `ref::RefParticle` – reference particle (energy at which all maps are built).
- `harmonic::Int` – h used to define the RF phase φ = −h ω0 Δt in laboratory coordinates.

Lattices are immutable; editing helpers ([`insert_cavity`](@ref),
[`with_reference`](@ref), ...) return a new Lattice.
"""
struct Lattice
    name::String
    elements::Vector{AbstractElement}
    ref::RefParticle
    harmonic::Int
end

"""
    path_length(lattice) -> Float64 [m]

Sum of element lengths: the length of the reference orbit. This is an arc
length, not necessarily the circumference of a closed ring.
"""
path_length(lattice::Lattice) = sum(elem_length, lattice.elements)
"""
    total_bend_angle(lattice) -> Float64 [rad]

Sum of the bending angles of all [`SBend`](@ref)s; 2π for a closed ring.
Used as a consistency check after loading/converting a lattice.
"""
total_bend_angle(lattice::Lattice) = sum(e -> e isa SBend ? e.angle : 0.0, lattice.elements)

"""
    unique_elements(lattice) -> Vector{AbstractElement}

Distinct elements (compared by value) in order of first appearance. Used to
print element tables and to build each map only once.
"""
function unique_elements(lattice::Lattice)
    seen = AbstractElement[]
    for e in lattice.elements
        any(isequal(e), seen) || push!(seen, e)
    end
    return seen
end

"""
    insert_cavity(lattice, cav::RFCavity; in_element) -> (Lattice, index)

Split the first drift named `in_element` into drift(L/2) – thin `cav` – drift(L/2),
keeping the path length. `in_element` is required. The cavity has zero length.
Returns the new lattice and the index of the cavity in it.
"""
function insert_cavity(lattice::Lattice, cav::RFCavity; in_element::AbstractString)
    i = findfirst(e -> e.name == in_element, lattice.elements)
    i === nothing && error("no element named $in_element")
    e = lattice.elements[i]
    e isa Drift || error("insert_cavity: $in_element is not a Drift")
    Ld = e.L / 2
    elems = vcat(lattice.elements[1:i-1], [Drift(e.name * "_u", Ld), cav, Drift(e.name * "_d", Ld)], lattice.elements[i+1:end])
    return Lattice(lattice.name * "+RF", elems, lattice.ref, cav.harmonic), i + 1
end

"""
    set_cavity_phase(lattice, phi_s) -> Lattice

Same lattice with every thin [`RFCavity`](@ref) set to synchronous phase `phi_s` [rad]
(used by `accelerate` to switch to a flat top).
"""
set_cavity_phase(lattice::Lattice, phi_s) =
    Lattice(lattice.name, [e isa RFCavity ? RFCavity(e.name, e.voltage_MV, e.harmonic, phi_s) : e for e in lattice.elements],
            lattice.ref, lattice.harmonic)

"""
    with_reference(lattice, ref) -> Lattice

Same elements (same normalised strengths) with another reference particle,
i.e. the magnets scaled with p0. Maps are rebuilt automatically (caches are keyed by energy).
"""
with_reference(lattice::Lattice, ref::RefParticle) = Lattice(lattice.name, lattice.elements, ref, lattice.harmonic)

# =============================================================================
# Structure-file parser
# =============================================================================

"""
    StructureFile(path, parameters, lattice)

Container for a parsed structure file: `path::String`, `parameters::Dict{String,Float64}`
(evaluated parameter section), `lattice::Lattice`.
"""
struct StructureFile
    path::String
    parameters::Dict{String,Float64}
    lattice::Lattice
end

_strip_comment(l) = (i = findfirst('#', l); i === nothing ? String(l) : String(l[1:prevind(l, i)]))

const _PARTICLES = Dict(
    "deuteron" => deuteron,
    "proton"   => proton,
    "electron" => electron,
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
        return SBend(name, L, get(kv, "angle", 0.0), get(kv, "k1", 0.0),
                     get(kv, "e1", 0.0), get(kv, "e2", 0.0), get(kv, "gap", 0.0), get(kv, "fint", 0.0))
    elseif k == "quadrupole" || k == "quad"
        return Quadrupole(name, L, get(kv, "k1", 0.0))
    elseif k == "sextupole" || k == "sext"
        return Sextupole(name, L, get(kv, "k2", 0.0))
    elseif k == "rfcavity" || k == "rf" || k == "cavity"
        V = haskey(kv, "v") ? kv["v"] : get(kv, "voltage_mv", 0.0)
        h = Int(round(haskey(kv, "h") ? kv["h"] : get(kv, "harmonic", 1.0)))
        ph = get(kv, "phi_s", get(kv, "phis", 0.0))
        return RFCavity(name, V, h, ph)
    elseif k == "marker"
        return Marker(name)
    else
        error("unknown element kind '$kind' for $name")
    end
end

function _expand_seq(tokens, elements, sequences; stack = String[])
    expanded = AbstractElement[]
    i = 1
    while i <= length(tokens)
        t = tokens[i]
        if i + 1 <= length(tokens) && tokens[i+1] == "*"
            n = Int(round(parse(Float64, t)))
            i + 2 <= length(tokens) || error("expected name after '*'")
            name = tokens[i+2]
            piece = _expand_name(name, elements, sequences, stack)
            for _ in 1:n
                append!(expanded, piece)
            end
            i += 3
        elseif t == "*"
            error("'*' without a count")
        else
            append!(expanded, _expand_name(t, elements, sequences, stack))
            i += 1
        end
    end
    return expanded
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

Read a Bell structure file (`.bell`).

# File format
Line-oriented; `#` starts a comment; keywords are case-insensitive.
```
reference
  particle     deuteron        # deuteron | proton | electron
  kinetic_MeV  270             # expression allowed
  harmonic     1               # optional (default 1); mass_MeV optional
parameters
  LB = 5.33
  ANG = 360/16 deg             # trailing 'deg' converts to rad
elements
  OBM drift      L = 2.0
  BM  sbend      L = LB  angle = ANG  k1 = 0  e1 = 0  e2 = 0  gap = 0  fint = 0
  QF  quadrupole L = 0.21 k1 = 0.5
  SD  sextupole  L = 0.1  k2 = 0
  RF  rfcavity   V = 0.01 h = 1 phi_s = 0
  CV  cavity     L = 0.5 V = 0.01 h = 1 phi_s = 0
  M1  marker
sequence cell                   # named sequence
  OBM BM QF ...
sequence                        # unnamed sequence = the ring
  8 * cell
```

# Algorithm
1. Each line is assigned to the current section.
2. Parameter lines `NAME = EXPR` are evaluated in order. Expressions allow
   `+ - * /`, parentheses, numbers and previously defined parameters.
3. Element cards are `NAME KIND key = expr ...`; keys are lower-cased, and
   missing keys default to 0. Kind aliases: `bend`, `quad`, `sext`, `rf`.
4. Sequences are expanded recursively: `N * name` repeats an element or a
   sequence, and cyclic references raise an error.
5. Build the `RefParticle` from particle/kinetic_MeV/mass_MeV and return a
   [`Lattice`](@ref) named after the file.

Errors are thrown (not warnings) for unknown kinds, unknown names, missing
reference data, or lines outside sections.
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
    return Lattice(name, elems, ref, harm)
end

_fmt(x::Real) = @sprintf("%.16g", x)

function _card_line(e::AbstractElement)
    nm = rpad(e.name, 8)
    if e isa Drift
        return @sprintf("%s %-12s L = %-20s", nm, "drift", _fmt(e.L))
    elseif e isa SBend
        s = @sprintf("%s %-12s L = %-20s angle = %-20s k1 = %s",
                     nm, "sbend", _fmt(e.L), _fmt(e.angle), _fmt(e.k1))
        for (k, v) in (("e1", e.e1), ("e2", e.e2), ("gap", e.gap), ("fint", e.fint))
            iszero(v) || (s *= "  $k = $(_fmt(v))")
        end
        return s
    elseif e isa Quadrupole
        return @sprintf("%s %-12s L = %-20s k1 = %s", nm, "quadrupole", _fmt(e.L), _fmt(e.k1))
    elseif e isa Sextupole
        return @sprintf("%s %-12s L = %-20s k2 = %s", nm, "sextupole", _fmt(e.L), _fmt(e.k2))
    elseif e isa RFCavity
        return @sprintf("%s %-12s V = %-20s h = %-6d phi_s = %s",
                        nm, "rfcavity", _fmt(e.voltage_MV), e.harmonic, _fmt(e.phi_s))
    elseif e isa Marker
        return @sprintf("%-8s marker", e.name)
    else
        error("cannot write element $(typeof(e))")
    end
end

"""
    write_structure(path, lattice; parameters = Dict(), extra_header = "") -> path

Write `lattice` in the format read by [`load_structure`](@ref): reference block,
optional `parameters`, one element card per distinct element (numbers with
`%.16g`, so a round trip is exact), and the ring sequence. `extra_header` is
copied as comment lines on top.
"""
function write_structure(path::AbstractString, lattice::Lattice;
                         parameters::Dict{String,Float64} = Dict{String,Float64}(),
                         extra_header = "")
    io = IOBuffer()
    println(io, "# Bell structure file")
    println(io, "# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>")
    isempty(extra_header) || println(io, extra_header)
    println(io)
    println(io, "reference")
    println(io, "  particle ", lattice.ref.particle.name)
    println(io, "  kinetic_MeV ", lattice.ref.Ekin_MeV)
    println(io, "  harmonic ", lattice.harmonic)
    if lattice.ref.mass_MeV != lattice.ref.particle.mass_MeV
        println(io, "  mass_MeV ", lattice.ref.mass_MeV)
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
    for e in lattice.elements
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
