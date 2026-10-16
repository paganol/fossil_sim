"""Spacecraft orbit: position and velocity around the Sun-Earth L2 point.

All vectors are barycentric, in the mean ecliptic frame of J2000
(astropy ``BarycentricMeanEcliptic(equinox="J2000")``, the frame of healpy's
``"E"`` coordinates). Times are seconds since ``mission_start`` (TDB).

The model follows the one of litebird_sim (``spacecraft.py``): L2 is the
Earth position and velocity scaled along the barycenter-Earth line, plus a
Lissajous orbit with an along-track and an out-of-ecliptic component.
Limitation: there is no radial component, while real L2 halo/Lissajous
orbits have one of comparable amplitude (~0.1 km/s in velocity, relevant
for the orbital dipole). A reference orbit from the mission team can
replace :class:`L2Orbit` by providing the same methods.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import (
    ICRS,
    BarycentricMeanEcliptic,
    CartesianRepresentation,
    get_body_barycentric,
    get_body_barycentric_posvel,
    solar_system_ephemeris,
)
from astropy.time import Time, TimeDelta
from scipy.interpolate import CubicSpline

from .config import OrbitConfig
from .constants import SECONDS_PER_YEAR


def icrs_to_ecliptic_matrix() -> np.ndarray:
    """Rotation matrix from ICRS to the J2000 mean ecliptic frame.

    Both frames are barycentric, so the transformation is a fixed rotation:
    it is obtained by transforming the three ICRS axes with astropy.

    >>> m = icrs_to_ecliptic_matrix()
    >>> bool(np.allclose(m @ m.T, np.eye(3)))
    True
    """
    axes = CartesianRepresentation(np.eye(3), unit=u.one)
    ecl = ICRS(axes).transform_to(BarycentricMeanEcliptic(equinox="J2000"))
    return ecl.cartesian.xyz.value  # columns are the images of the ICRS axes


class L2Orbit:
    """Approximate L2 orbit with an optional Lissajous term.

    Ephemerides are evaluated with astropy on a regular grid covering
    ``[0, duration_s]`` and interpolated with cubic splines.

    Parameters
    ----------
    config : OrbitConfig
        Orbit parameters.
    mission_start : astropy.time.Time
        Time corresponding to ``time_s = 0``.
    duration_s : float
        Length of the time span to cover, in seconds.
    grid_step_s : float
        Spacing of the ephemeris grid, in seconds.
    """

    def __init__(
        self,
        config: OrbitConfig,
        mission_start: Time,
        duration_s: float,
        grid_step_s: float = 3600.0,
    ) -> None:
        self.config = config
        num_points = int(np.ceil(duration_s / grid_step_s)) + 4
        grid_s = (np.arange(num_points) - 1) * grid_step_s  # margin on both sides
        times = mission_start + TimeDelta(grid_s, format="sec")
        rot = icrs_to_ecliptic_matrix()
        with solar_system_ephemeris.set("builtin"):
            earth_pos, earth_vel = get_body_barycentric_posvel("earth", times)
            sun_pos = get_body_barycentric("sun", times)
        earth_pos_km = rot @ earth_pos.xyz.to_value(u.km)
        earth_vel_km_s = rot @ earth_vel.xyz.to_value(u.km / u.s)
        sun_pos_km = rot @ sun_pos.xyz.to_value(u.km)
        self._earth_pos = CubicSpline(grid_s, earth_pos_km.T)
        self._earth_vel = CubicSpline(grid_s, earth_vel_km_s.T)
        self._sun_pos = CubicSpline(grid_s, sun_pos_km.T)

    def position_velocity(self, time_s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Spacecraft position and velocity.

        Parameters
        ----------
        time_s : array_like
            Times in seconds since ``mission_start``.

        Returns
        -------
        pos_km, vel_km_s : ndarray, shape (n, 3)
            Barycentric position (km) and velocity (km/s), J2000 ecliptic.
        """
        time_s = np.atleast_1d(np.asarray(time_s, dtype=np.float64))
        earth_pos_km = self._earth_pos(time_s)
        earth_vel_km_s = self._earth_vel(time_s)
        earth_dist_km = np.linalg.norm(earth_pos_km, axis=1, keepdims=True)
        scale = 1.0 + self.config.earth_l2_distance_km / earth_dist_km
        pos_km = earth_pos_km * scale
        # d/dt of (r + d_L2 r̂): the radial part of v is unchanged, the tangential
        # part is scaled; exact derivative of pos_km above
        radial_dir = earth_pos_km / earth_dist_km
        v_radial = np.sum(earth_vel_km_s * radial_dir, axis=1, keepdims=True)
        vel_km_s = earth_vel_km_s * scale + v_radial * radial_dir * (1.0 - scale)

        lissajous_pos_km, lissajous_vel_km_s = self._lissajous(
            time_s, earth_pos_km, earth_vel_km_s
        )
        return pos_km + lissajous_pos_km, vel_km_s + lissajous_vel_km_s

    def _lissajous(self, time_s, earth_pos_km, earth_vel_km_s):
        """Lissajous displacement: along-track (radius 1) and ecliptic-z (radius 2)."""
        cfg = self.config
        omega1_rad_s = 2 * np.pi * cfg.lissajous_freq1_cycles_yr / SECONDS_PER_YEAR
        omega2_rad_s = 2 * np.pi * cfg.lissajous_freq2_cycles_yr / SECONDS_PER_YEAR
        phase1_rad = omega1_rad_s * time_s[:, None]
        phase2_rad = omega2_rad_s * time_s[:, None] + np.deg2rad(
            cfg.lissajous_phase_deg
        )

        # Along-track unit vector t̂ = ẑ × r̂ and its time derivative
        earth_dist_km = np.linalg.norm(earth_pos_km, axis=1, keepdims=True)
        radial_dir = earth_pos_km / earth_dist_km
        radial_dir_dot = (
            earth_vel_km_s
            - np.sum(earth_vel_km_s * radial_dir, axis=1, keepdims=True) * radial_dir
        ) / earth_dist_km
        z_axis = np.array([0.0, 0.0, 1.0])
        track_dir = np.cross(z_axis, radial_dir)
        track_dir_dot = np.cross(z_axis, radial_dir_dot)

        r1, r2 = cfg.lissajous_radius1_km, cfg.lissajous_radius2_km
        pos_km = r1 * np.cos(phase1_rad) * track_dir + r2 * np.sin(phase2_rad) * z_axis
        vel_km_s = (
            -r1 * omega1_rad_s * np.sin(phase1_rad) * track_dir
            + r1 * np.cos(phase1_rad) * track_dir_dot
            + r2 * omega2_rad_s * np.cos(phase2_rad) * z_axis
        )
        return pos_km, vel_km_s

    def anti_sun_direction(self, time_s: np.ndarray) -> np.ndarray:
        """Unit vector from the Sun to the spacecraft, shape (n, 3)."""
        time_s = np.atleast_1d(np.asarray(time_s, dtype=np.float64))
        pos_km, _ = self.position_velocity(time_s)
        direction = pos_km - self._sun_pos(time_s)
        return direction / np.linalg.norm(direction, axis=1, keepdims=True)

    def anti_sun_longitude_rad(self, time_s: np.ndarray) -> np.ndarray:
        """Ecliptic longitude of the anti-Sun direction, unwrapped, in radians."""
        direction = self.anti_sun_direction(time_s)
        return np.unwrap(np.arctan2(direction[:, 1], direction[:, 0]))
