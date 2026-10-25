"""FOSSIL scanning strategy: Spin axis per Ring and boresight per IFG.

The geometry is the one of litebird_sim's ``SpinningScanningStrategy``
(spin, cone around the anti-Sun direction, yearly motion), with one
difference: the Spin-axis longitude and the Precession phase are frozen
within each Ring, so that its Scans repeat the same circle.

Conventions:

- Spin-axis longitude from the Ring stepping rule (CONTEXT.md):

  * ``equal_per_cycle``: per Cooler cycle the axis advances by the actual
    solar motion over the cycle, in ``rings_per_cycle`` equal steps, with an
    offset that centres the Solar Aspect Angle excursion;
  * ``anchor_ring_start``: anti-Sun longitude at the start of each Ring;
  * ``fixed``: anti-Sun longitude at mission start plus ``fixed_ring_step_deg``
    per Ring (no re-anchoring: test/diagnostic only).

  During the Recycling window the axis stays at the last Ring's position.
- Precession: the axis sits on a cone of opening ``precession_amplitude_deg``
  around the nominal (ecliptic-plane) axis, with phase
  ``2π t_ring_mid / P + phase0`` (0 = towards ecliptic east, π/2 = north).
- Spin phase 0 (start of every Scan) points at the point of the circle
  closest to the north ecliptic pole; the spin is right-handed about the axis.
"""

from dataclasses import dataclass

import numpy as np
from numba import njit

from .config import Config
from .operations import CycleTimeline
from .orbit import L2Orbit


@dataclass
class CyclePointing:
    """Pointing of one Cooler cycle."""

    spin_axis: np.ndarray  # (rings_per_cycle, 3), unit vectors
    boresight: np.ndarray  # (n_ifg, 3), unit vectors at the IFG mid-time
    saa_rad: np.ndarray  # (n_ifg,), Spin axis / anti-Sun angle


def equal_step_offset(
    scans_per_ring: int, rings_per_cycle: int, scans_per_cycle: int
) -> float:
    """Offset ``o`` (in units of steps) of the ``equal_per_cycle`` rule.

    Ring ``r`` sits at ``λ0 + (r + o) Δ / n``, where ``Δ`` is the solar motion
    over the cycle. With a uniform solar motion, ``o`` centres the anti-Sun
    excursion around the axis over the whole cycle (Recycling included).

    >>> round(equal_step_offset(23, 3, 87), 4)
    0.2931
    """
    n = rings_per_cycle
    f = scans_per_ring / scans_per_cycle  # fraction of the cycle per Ring
    return 0.5 * ((n - 1) * (n * f - 1) + 1)


def ring_axis_longitudes_rad(
    config: Config, orbit: L2Orbit, cycle_id: int
) -> np.ndarray:
    """Nominal Spin-axis ecliptic longitude of each Ring of a cycle."""
    tl, sc = config.timeline, config.scanning
    cycle_start_s = cycle_id * tl.cycle_duration_s
    ring_duration_s = tl.scans_per_ring * tl.scan_duration_s
    rings = np.arange(tl.rings_per_cycle)
    if sc.ring_stepping == "equal_per_cycle":
        start_end = orbit.anti_sun_longitude_rad(
            np.array([cycle_start_s, cycle_start_s + tl.cycle_duration_s])
        )
        step_rad = (start_end[1] - start_end[0]) / tl.rings_per_cycle
        offset = equal_step_offset(
            tl.scans_per_ring, tl.rings_per_cycle, tl.scans_per_cycle
        )
        return start_end[0] + (rings + offset) * step_rad
    if sc.ring_stepping == "anchor_ring_start":
        return orbit.anti_sun_longitude_rad(cycle_start_s + rings * ring_duration_s)
    # "fixed"
    first_rad = orbit.anti_sun_longitude_rad(np.array([0.0]))[0]
    ring_id = cycle_id * tl.rings_per_cycle + rings
    return first_rad + ring_id * np.deg2rad(sc.fixed_ring_step_deg)


def ring_spin_axes(config: Config, orbit: L2Orbit, cycle_id: int) -> np.ndarray:
    """Spin axis (unit vector, J2000 ecliptic) of each Ring of a cycle."""
    tl, sc = config.timeline, config.scanning
    longitude_rad = ring_axis_longitudes_rad(config, orbit, cycle_id)
    nominal = np.stack(
        [np.cos(longitude_rad), np.sin(longitude_rad), np.zeros_like(longitude_rad)],
        axis=1,
    )
    amplitude_rad = np.deg2rad(sc.precession_amplitude_deg)
    if amplitude_rad == 0.0:
        return nominal
    ring_duration_s = tl.scans_per_ring * tl.scan_duration_s
    ring_mid_s = (
        cycle_id * tl.cycle_duration_s
        + (np.arange(tl.rings_per_cycle) + 0.5) * ring_duration_s
    )
    phase_rad = 2 * np.pi * ring_mid_s / (
        sc.precession_period_days * 86400.0
    ) + np.deg2rad(sc.precession_phase_deg)
    east = np.stack(
        [-np.sin(longitude_rad), np.cos(longitude_rad), np.zeros_like(longitude_rad)],
        axis=1,
    )
    north = np.array([0.0, 0.0, 1.0])
    offset_dir = np.cos(phase_rad)[:, None] * east + np.sin(phase_rad)[:, None] * north
    return np.cos(amplitude_rad) * nominal + np.sin(amplitude_rad) * offset_dir


@njit(cache=True)
def _boresight_kernel(axes, ring_of_ifg, spin_phase_rad, spin_boresight_angle_rad, out):
    """Boresight per IFG; phase 0 is the circle point closest to the ecliptic north pole."""
    cos_b = np.cos(spin_boresight_angle_rad)
    sin_b = np.sin(spin_boresight_angle_rad)
    for i in range(spin_phase_rad.size):
        ax, ay, az = axes[ring_of_ifg[i]]
        # u: north pole projected orthogonally to the axis; w = axis × u
        ux, uy, uz = -az * ax, -az * ay, 1.0 - az * az
        norm = np.sqrt(ux * ux + uy * uy + uz * uz)
        ux, uy, uz = ux / norm, uy / norm, uz / norm
        wx, wy, wz = ay * uz - az * uy, az * ux - ax * uz, ax * uy - ay * ux
        c, s = np.cos(spin_phase_rad[i]), np.sin(spin_phase_rad[i])
        out[i, 0] = cos_b * ax + sin_b * (c * ux + s * wx)
        out[i, 1] = cos_b * ay + sin_b * (c * uy + s * wy)
        out[i, 2] = cos_b * az + sin_b * (c * uz + s * wz)


def cycle_pointing(
    config: Config, orbit: L2Orbit, timeline: CycleTimeline
) -> CyclePointing:
    """Spin axes, boresight and Solar Aspect Angle for one Cooler cycle."""
    axes = ring_spin_axes(config, orbit, timeline.cycle_id)
    ring = timeline.ring_in_cycle.astype(np.int64)
    spin_phase_rad = (
        2 * np.pi * (timeline.ifg_in_scan + 0.5) / config.timeline.ifgs_per_scan
    )
    boresight = np.empty((ring.size, 3))
    _boresight_kernel(
        axes,
        ring,
        spin_phase_rad,
        np.deg2rad(config.scanning.spin_boresight_angle_deg),
        boresight,
    )
    anti_sun = orbit.anti_sun_direction(timeline.time_mid_s)
    cos_saa = np.clip(np.sum(axes[ring] * anti_sun, axis=1), -1.0, 1.0)
    return CyclePointing(
        spin_axis=axes, boresight=boresight, saa_rad=np.arccos(cos_saa)
    )


def ring_circles(
    config: Config,
    axes: np.ndarray,
    points_per_ifg: int = 1,
    phase_offset: float = 0.5,
) -> np.ndarray:
    """Boresight along the circle of each Ring, shape (n_rings, n_points, 3).

    The Scans of a Ring retrace this circle (Spin axis frozen within the Ring).
    Point g is at spin phase 2π (g + ``phase_offset``) / n_points, with
    n_points = ``points_per_ifg`` × ifgs_per_scan; the defaults give the IFG
    centres.
    """
    n_ifg = config.timeline.ifgs_per_scan * points_per_ifg
    phase_rad = 2 * np.pi * (np.arange(n_ifg) + phase_offset) / n_ifg
    ring = np.repeat(np.arange(axes.shape[0]), n_ifg)
    out = np.empty((ring.size, 3))
    _boresight_kernel(
        axes,
        ring,
        np.tile(phase_rad, axes.shape[0]),
        np.deg2rad(config.scanning.spin_boresight_angle_deg),
        out,
    )
    return out.reshape(axes.shape[0], n_ifg, 3)
