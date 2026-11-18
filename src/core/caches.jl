# Registry for map memoisation. Included into Bell with core.
# The dictionaries live in the map modules. Bell.__init__ registers an
# emptier for each one after those modules exist. This file does not name them.

const _CACHE_EMPTIERS = Vector{Function}()

"""
    register_cache_empty!(empty!)

Record a zero-argument function that [`clear_caches!`](@ref) will call.
`Bell` registers the Lie, DA, linear, and fringe caches when it loads.
"""
function register_cache_empty!(empty!::Function)
    push!(_CACHE_EMPTIERS, empty!)
    return nothing
end

"""
    clear_caches!()

Empty the memoised maps and the shared DA spaces registered at load:

- DA element maps (`DATracking._DACACHE`) and shared DA spaces (`_SPACES`);
- Lie element maps (`Tracking._MCACHE`) and Taylor element maps (`_TMCACHE`);
- linear `SA` matrices (`Linear._HCACHE`);
- fringe edge maps (`FringeMaps._ECACHE`).

The next map, matrix or edge rebuilds the entry it needs. Objects already
returned stay as they are; they are not updated in place.
"""
function clear_caches!()
    for empty! in _CACHE_EMPTIERS
        empty!()
    end
    return nothing
end
