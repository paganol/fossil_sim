# Conventions

## Units

Measure units are always explicit in variable, parameter and configuration-key
names, as a lowercase suffix: `freq_ghz`, `time_s`, `temperature_k`,
`theta_rad`, `v_km_s`, `intensity_mjy_sr`. Internally fossil_sim uses

| Quantity | Unit |
|---|---|
| frequency | GHz |
| intensity | MJy sr⁻¹ (1 MJy sr⁻¹ = 10⁻²⁰ W m⁻² Hz⁻¹ sr⁻¹) |
| time | s since `mission_start` |
| angles | rad (degrees only in the configuration, suffix `_deg`) |
| temperature | K |
| distance, velocity | km, km s⁻¹ |

## Time

Absolute times are represented as seconds since `mission_start`, a TDB epoch
given in the configuration (the start of science operations at L2). A float64
second count keeps a resolution of about 20 ns over the 3.5-year mission,
whereas astropy `Time` objects are used only to evaluate ephemerides.

## Reference frame

All vectors (pointing, positions, velocities) are expressed in the barycentric
mean ecliptic frame of J2000, astropy `BarycentricMeanEcliptic(equinox="J2000")`.
This is the frame of healpy's `"E"` coordinates, so maps produced by
fossil_sim and sky models rotated with `healpy.Rotator(coord=[..., "E"])`
agree without further rotation. Maps use HEALPix RING ordering.

## Reproducibility

The configuration holds a single master seed. Random streams are derived from
it by content, `SeedSequence(master, spawn_key=(component, cycle_id, du))`,
never by MPI rank: a run gives identical results on any number of processes,
and any Cooler cycle can be regenerated in isolation
({doc}`../adr/0001-stateless-cooler-cycles`).

## Flags

Each IFG carries an 8-bit flag word; only IFGs with no bit set (Valid IFGs)
enter the maps.

| Bit | Value | Meaning |
|---|---|---|
| `FLAG_RECYCLING` | 1 | inside a Recycling window |
| `FLAG_MODE_SWITCH` | 2 | dead time after a Mode change |
| `FLAG_SAA_EXCEEDED` | 4 | Solar Aspect Angle above `max_saa_deg` |
| `FLAG_SLEW` | 8 | reserved (slews are instantaneous) |
