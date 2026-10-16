# Cooler cycles are independent, stateless units

The mission is simulated one Cooler cycle at a time, and no state is carried from one cycle to the next: every quantity (pointing, Mode, BBIR set-point, noise) is a deterministic function of global indices and absolute time. Random streams are keyed by content, `SeedSequence(master, spawn_key=(component, cycle_id, du))`, not by MPI rank. This makes cycles embarrassingly parallel (rank r takes cycles r, r+N, …, maps are summed at the end), keeps results identical for any number of ranks, and lets any single cycle be regenerated in isolation.

## Considered Options

- Per-rank RNG hierarchy (as in litebird_sim): rejected, because results then depend on the number of ranks.
- Carrying state across cycles (continuous 1/f, thermal drifts): rejected for now.

## Consequences

Processes with memory longer than a cycle (1/f noise, BBIR thermal transients, slow drifts) must either reset during the Recycling window, where data are flagged anyway, or be modelled as deterministic functions of absolute time. Revisit this ADR if a systematic needs true cross-cycle memory.
