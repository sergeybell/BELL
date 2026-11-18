# tracking

Runtime state and dispatch over the map engines in `src/maps/`. `state.jl` (`BeamState`) is `ParticleState` and `Beam`: canonical coordinates and a spin slot that is only stored. `trackers.jl` is `FieldTracker`, `LieTracker`, and `MapTracker`. Each one calls an existing function (`exact_track_turn`, `track` on factors, or a one-turn map).

`BeamState` imports no Bell module. `Trackers` depends on `BeamState`, `Structure`, `Integrators`, `FieldTrack`, Lie `Tracking`, and `LieAlgebra`. The `Bell` facade re-exports both. Diagram: [src/README.md](../README.md).

A new tracker that only selects an existing engine goes in `trackers.jl`. New particle data goes in `state.jl`. A new equation goes in `src/physics/` or `src/maps/`.
