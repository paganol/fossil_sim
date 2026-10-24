"""Spectral maps for the first months of the mission, with a markdown report.

Runs the baseline twice with the same seed, with and without noise: the
difference is the noise realization alone, which must match the variance maps.

Usage::

    uv run python examples/spectral_maps.py [config.toml] [output_dir] [days]
"""

import sys
import time
from dataclasses import replace
from pathlib import Path

import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fossil_sim.config import load_config
from fossil_sim.simulation import run_simulation

REPO = Path(__file__).parents[1]
SHOWN_FREQ_GHZ = (105.0, 405.0, 1005.0)


def main(config_path: Path, output_dir: Path, days: float) -> None:
    config = load_config(config_path)
    tl = config.timeline
    num_of_cycles = int(np.ceil(days * 86400 / tl.cycle_duration_s))
    output_dir.mkdir(parents=True, exist_ok=True)

    start = time.perf_counter()
    noisy = run_simulation(config, num_of_cycles=num_of_cycles)
    quiet = run_simulation(
        replace(config, noise=replace(config.noise, enabled=False)),
        num_of_cycles=num_of_cycles,
    )
    elapsed_s = time.perf_counter() - start
    noisy.spectral_maps.save(output_dir / "maps")
    sm, freq_ghz = noisy.spectral_maps, noisy.spectral_maps.freq_ghz

    lost = noisy.lost_fraction()
    lines = [
        "# FOSSIL spectral maps",
        "",
        (
            f"Config: `{config_path}` (hash `{config.hash()[:12]}`, "
            f"git `{config.provenance()['git_sha'][:10]}`)"
        ),
        "",
        (
            f"- {num_of_cycles} Cooler cycles ({num_of_cycles * tl.cycle_duration_s / 86400:.1f} d), "
            f"two runs (noisy + noiseless) in {elapsed_s / 60:.1f} min"
        ),
        f"- IFGs: {noisy.num_ifgs:,}, valid: {noisy.num_valid_ifgs:,}",
        "- Time lost: "
        + ", ".join(f"{name} {100 * frac:.2f}%" for name, frac in lost.items()),
        f"- Max Solar Aspect Angle: {noisy.max_saa_deg:.3f} deg",
        f"- Frequency grid: {freq_ghz[0]:.0f}–{freq_ghz[-1]:.0f} GHz, {freq_ghz.size} bins",
        (
            f"- Static sky: {', '.join(config.static_sky.models) or 'none'}; "
            f"zodi: {config.zodi.model if config.zodi.enabled else 'off'}"
        ),
        "",
        "## Noise validation",
        "",
        (
            "χ = (noisy − noiseless) / √variance over the observed pixels; "
            "its standard deviation must be 1 in every bin."
        ),
        "",
        "| nside | sky coverage | hits per pixel (median) | std(χ) min | std(χ) max |",
        "|---|---|---|---|---|",
    ]
    for nside in config.maps.nsides:
        seen = sm.sky_count[nside] > 0
        variance = sm.sky_variance_map(nside)[:, seen]
        residual = (
            sm.sky_map(nside)[:, seen] - quiet.spectral_maps.sky_map(nside)[:, seen]
        )
        chi_std = np.std(residual / np.sqrt(variance), axis=1)
        lines.append(
            f"| {nside} | {100 * seen.mean():.2f}% | {np.median(sm.sky_count[nside][seen]):.0f} "
            f"| {chi_std.min():.3f} | {chi_std.max():.3f} |"
        )

    nside = min(config.maps.nsides)
    seen = sm.sky_count[nside] > 0
    sigma = np.sqrt(sm.sky_variance_map(nside)[:, seen])
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.semilogy(freq_ghz, np.median(sigma, axis=1) * 1e3)
    ax.set(
        xlabel="GHz",
        ylabel="median σ per pixel [kJy/sr]",
        title=f"White noise, nside {nside}, {days:.0f} d",
    )
    fig.tight_layout()
    fig.savefig(output_dir / "noise_spectrum.png", dpi=90)
    plt.close(fig)
    lines += ["", "![](noise_spectrum.png)", "", f"## Maps (nside {nside})", ""]

    sky = sm.sky_map(nside)
    std_map = np.where(seen, np.sqrt(np.abs(sm.sky_variance_map(nside))), hp.UNSEEN)
    for f in SHOWN_FREQ_GHZ:
        i = int(np.argmin(np.abs(freq_ghz - f)))
        fig = plt.figure(figsize=(14, 4))
        hp.mollview(
            np.where(seen, sky[i] - np.median(sky[i][seen]), hp.UNSEEN),
            fig=fig,
            sub=(1, 3, 1),
            unit="MJy/sr",
            norm="hist",
            title=f"sky − median, {freq_ghz[i]:.0f} GHz",
        )
        hp.mollview(
            np.where(seen, sky[i] - quiet.spectral_maps.sky_map(nside)[i], hp.UNSEEN),
            fig=fig,
            sub=(1, 3, 2),
            unit="MJy/sr",
            title="noise (noisy − noiseless)",
        )
        hp.mollview(
            std_map[i],
            fig=fig,
            sub=(1, 3, 3),
            unit="MJy/sr",
            title="σ per pixel",
        )
        png = output_dir / f"maps_{freq_ghz[i]:04.0f}GHz.png"
        fig.savefig(png, dpi=90)
        plt.close(fig)
        lines.append(f"![]({png.name})")

    report = output_dir / "report.md"
    report.write_text("\n".join(lines) + "\n")
    print(report.read_text())


if __name__ == "__main__":
    config_path = (
        Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "configs" / "baseline.toml"
    )
    output_dir = (
        Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "output" / "spectral_maps"
    )
    days = float(sys.argv[3]) if len(sys.argv) > 3 else 183.0
    main(config_path, output_dir, days)
