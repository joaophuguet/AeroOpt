"""Executable geometry-only demonstration. Use a fresh output directory."""

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aeroopt.geometry import MultiElementAirfoil
from aeroopt.geometry.plotting import compare_configurations, save_figure


def main() -> None:
    """Import, normalize, place, reject/export, plot, and summarize candidates."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "demo")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output directory must be empty: {output}. Choose --output with a new path.")
    base = MultiElementAirfoil.from_dat(
        ROOT / "examples/data/demo_main.dat", ROOT / "examples/data/demo_flap.dat",
        flap_chord_ratio=0.30, flap_reference_position=(0.95, -0.02),
    )
    # Intentional valid, intersecting, and low-clearance cases; no performance claims.
    designs = [("downstream", 0.0, 0.15, 0.0),
               ("deflected", 25.0, 0.08, -0.04),
               ("lowered", 40.0, 0.0, -0.08),
               ("intersecting", 0.0, -0.20, 0.02),
               ("low_clearance", 0.0, 0.0505, 0.02),
               ("displaced", 15.0, 0.12, 0.04)]
    clearance = 0.002
    geometries, labels, summary = [], [], []
    output.mkdir(parents=True, exist_ok=True)
    for label, angle, dx, dy in designs:
        geometry = base.configure(flap_deflection_deg=angle, flap_dx=dx, flap_dy=dy)
        report = geometry.validate(min_clearance=clearance)
        geometries.append(geometry)
        labels.append(label)
        summary.append({"id": label, "parameters": geometry.parameters.to_dict(),
                        "validation": report.to_dict()})
        if report.valid:
            geometry.export(output / label, configuration_id=label, min_clearance=clearance)
        codes = ", ".join(issue.code for issue in report.issues) or "ok"
        distance = "undefined" if report.min_distance is None else f"{report.min_distance:.8f} c"
        print(f"{label:16} {'VALID' if report.valid else 'REJECTED':8} distance={distance}: {codes}")
    fig = compare_configurations(geometries, labels=labels, min_clearance=clearance)
    save_figure(fig, output / "comparison")
    plt.close(fig)
    fig = geometries[1].plot().figure
    save_figure(fig, output / "deflected_detail")
    plt.close(fig)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    valid = sum(item["validation"]["valid"] for item in summary)
    print(f"\n{valid}/{len(summary)} valid; {len(summary) - valid} rejected. Outputs: {output}")
    print("Synthetic demonstration contours only. No aerodynamic results are computed.")


if __name__ == "__main__":
    main()
