"""Hit maps for 30 days and 6 months, with a short markdown report.

Usage::

    uv run python examples/hitmaps.py [config.toml] [output_dir]
"""

import sys
import time
from pathlib import Path

import healpy as hp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fossil_sim.config import load_config
from fossil_sim.simulation import run_simulation

REPO = Path(__file__).parents[1]
CASES_DAYS = {"30d": 30.0, "6mo": 183.0}


def main(config_path: Path, output_dir: Path) -> None:
    config = load_config(config_path)
    tl = config.timeline
    lines = [
        "# FOSSIL hit maps",
        "",
        (
            f"Config: `{config_path}` (hash `{config.hash()[:12]}`, "
            f"git `{config.provenance()['git_sha'][:10]}`)"
        ),
        "",
        "## Timeline",
        "",
        (
            f"- Scan: {tl.scan_duration_s:.0f} s, {tl.ifgs_per_scan} IFGs, "
            f"spin period {tl.scan_duration_s / 60:.0f} min "
            f"({3600 / tl.scan_duration_s:.2f} rph; Step2 text says ~2 rph)"
        ),
        (
            f"- Step per IFG: {360 / tl.ifgs_per_scan:.2f} deg, "
            f"{config.instrument.beam_fwhm_deg / (360 / tl.ifgs_per_scan):.2f} samples per beam FWHM"
        ),
        f"- Ring: {tl.scans_per_ring} Scans = {tl.scans_per_ring * tl.scan_duration_s / 3600:.2f} h",
        (
            f"- Cooler cycle: {tl.rings_per_cycle} Rings + {tl.recycling_scans} Recycling Scans "
            f"= {tl.cycle_duration_s / 3600:.2f} h"
        ),
        "",
    ]
    for label, days in CASES_DAYS.items():
        num_of_cycles = int(np.ceil(days * 86400 / tl.cycle_duration_s))
        start = time.perf_counter()
        result = run_simulation(config, num_of_cycles=num_of_cycles)
        elapsed_s = time.perf_counter() - start
        case_dir = output_dir / label
        result.hit_maps.save(case_dir)
        lost = result.lost_fraction()
        lines += [
            f"## {label}: {num_of_cycles} Cooler cycles ({num_of_cycles * tl.cycle_duration_s / 86400:.1f} d)",
            "",
            f"- IFGs: {result.num_ifgs:,}, valid: {result.num_valid_ifgs:,} (run time {elapsed_s:.1f} s)",
            "- Time lost: "
            + ", ".join(f"{name} {100 * frac:.2f}%" for name, frac in lost.items()),
            "- Valid IFGs per Mode: "
            + ", ".join(f"{m} {n:,}" for m, n in result.mode_counts.items() if n),
            (
                f"- Max Solar Aspect Angle: {result.max_saa_deg:.3f} deg "
                f"(limit {config.scanning.max_saa_deg} deg)"
            ),
        ]
        for nside in config.maps.nsides:
            total = result.hit_maps.total(nside)
            coverage = np.mean(total > 0)
            first_s = result.hit_maps.first_hit_time_s[nside]
            full_day = np.max(first_s) / 86400 if np.isfinite(first_s).all() else None
            lines.append(
                f"- nside {nside}: sky coverage {100 * coverage:.2f}%"
                + (f", full sky reached at day {full_day:.1f}" if full_day else "")
            )
            fig = plt.figure(figsize=(14, 4))
            for i, (title, hits) in enumerate(
                [
                    ("total", total),
                    ("B10", result.hit_maps.hits[nside][2]),
                    ("B01", result.hit_maps.hits[nside][1]),
                ]
            ):
                hp.mollview(
                    hits,
                    fig=fig,
                    sub=(1, 3, i + 1),
                    title=f"{label} {title}, nside {nside}",
                    unit="hits",
                    norm="log" if hits.max() > 0 else None,
                    min=1,
                )
            png = case_dir / f"hits_nside{nside:04d}.png"
            fig.savefig(png, dpi=90)
            plt.close(fig)
            lines.append(f"  ![]({png.relative_to(output_dir)})")
        days_axis = np.linspace(0, num_of_cycles * tl.cycle_duration_s / 86400, 200)
        fig, ax = plt.subplots(figsize=(6, 3.5))
        for nside in config.maps.nsides:
            ax.plot(
                days_axis,
                result.hit_maps.coverage_fraction(nside, days_axis * 86400),
                label=f"nside {nside}",
            )
        ax.set(xlabel="days", ylabel="sky fraction with ≥1 valid hit", title=label)
        ax.legend()
        fig.tight_layout()
        png = case_dir / "coverage.png"
        fig.savefig(png, dpi=90)
        plt.close(fig)
        lines += [f"  ![]({png.relative_to(output_dir)})", ""]

    report = output_dir / "report.md"
    report.write_text("\n".join(lines) + "\n")
    print(report.read_text())


if __name__ == "__main__":
    config_path = (
        Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "configs" / "baseline.toml"
    )
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "output" / "hitmaps"
    main(config_path, output_dir)
