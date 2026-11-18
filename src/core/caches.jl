# Map memoisation shared by the Lie, DA and fringe paths.
# Included into Bell after those modules exist.

"""
    clear_caches!()

Empty the memoised maps and the shared DA spaces:

- DA element maps (`DATracking._DACACHE`) and shared DA spaces (`_SPACES`);
- Lie element maps (`Tracking._MCACHE`) and Taylor element maps (`_TMCACHE`);
- linear `SA` matrices (`Linear._HCACHE`);
- fringe edge maps (`FringeMaps._ECACHE`).

The next map, matrix or edge rebuilds the entry it needs. Objects already
returned stay as they are; they are not updated in place.
"""
function clear_caches!()
    empty!(DATracking._DACACHE)
    empty!(DATracking._SPACES)
    empty!(Tracking._MCACHE)
    empty!(Tracking._TMCACHE)
    empty!(Linear._HCACHE)
    empty!(FringeMaps._ECACHE)
    return nothing
end
