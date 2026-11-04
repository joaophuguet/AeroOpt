"""Geometry-only v0.1.1 validation; offline analytic references + optional UIUC.

Use --fetch-uiuc for opt-in downloads, or --reference-dir for preserved originals.
Downloaded coordinates have no redistribution permission assumed by this study.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapely

from aeroopt import __version__
from aeroopt.geometry import (Airfoil, FlapParameters, MultiElementAirfoil, naca4_points,
                              read_dat, validate_elements)
from aeroopt.geometry.plotting import compare_configurations, save_figure
from _naca_reference import continuous_metrics, continuous_slot_distance
from _reference_data import published_sources
from _validation_cases import validation_cases


def _json(path: Path, data: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def _table(output: Path, stem: str, rows: list[dict]) -> None:
    _json(output / f"{stem}.json", rows)
    with (output / f"{stem}.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _perimeter(points: np.ndarray) -> float:
    return float(np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1).sum())


def _signed_area(points: np.ndarray) -> float:
    p = points - points[0]
    return float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)) / 2)


def _invariance(foil: Airfoil) -> list[dict]:
    results = []
    before = foil.original_points.copy()
    base = foil.normalized_points
    distances = np.linalg.norm(base[:, None] - base[None, :], axis=2)
    for scale in (1e-6, 1., 1e6):
        for angle in (-75., 0., 138.):
            a = np.deg2rad(angle)
            rotation = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
            transformed = scale * (before @ rotation.T + [2.3, -1.7])
            moved = Airfoil.from_points(transformed, leading_edge_index=foil.leading_edge_index)
            points = moved.normalized_points
            distance_error = float(np.max(np.abs(
                np.linalg.norm(points[:, None] - points[None, :], axis=2) - distances)))
            point_error = float(np.max(np.linalg.norm(points - base, axis=1)))
            area_error = abs(moved.area / (foil.area * scale**2) - 1)
            orientation = np.sign(_signed_area(points)) == np.sign(_signed_area(base))
            if max(point_error, distance_error, area_error) > 2e-12 or not orientation:
                raise AssertionError("normalization invariance exceeded explicit 2e-12 study tolerance")
            np.testing.assert_array_equal(moved.original_points, transformed)
            results.append({"profile": foil.name, "scale": scale, "rotation_deg": angle,
                            "max_normalized_point_error": point_error,
                            "max_relative_distance_error": distance_error,
                            "relative_area_error": area_error,
                            "physical_orientation_preserved": bool(orientation)})
    np.testing.assert_array_equal(foil.original_points, before)
    return results


def run(output: Path, reference_dir: Path | None, fetch: bool) -> None:
    """Execute all diagnostics, fail on unexpected classification or round-trip error."""
    output = output.resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError("results must stay inside the AeroOpt workspace")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty {output}; choose a new --output")
    output.mkdir(parents=True, exist_ok=True)
    provenance, profiles, normalization = [], [], []
    generated: dict[tuple, Airfoil] = {}

    def analytic(code: str, count: int, closed: bool = False, spacing: str = "cosine") -> Airfoil:
        key = (code, count, closed, spacing)
        if key not in generated:
            path = output / "analytic" / f"naca{code}_{count}_{spacing}_{'closed' if closed else 'finite'}_te.dat"
            path.parent.mkdir(parents=True, exist_ok=True)
            points = naca4_points(code, count, closed_te=closed, spacing=spacing)
            with path.open("x", encoding="utf-8", newline="\n") as stream:
                np.savetxt(stream, points, fmt="%.17g",
                           header=f"AeroOpt_COORDINATES_v0.1 ANALYTIC_NACA_{code}", comments="")
            foil = Airfoil.from_dat(path, leading_edge_index=count // 2)
            generated[key] = foil
            provenance.append({"name": f"NACA {code}", "file": str(path.relative_to(output)),
                               "origin": "analytically generated geometry; not experimental data",
                               "point_count": count, "spacing": spacing, "closed_te": closed,
                               "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                               "coordinate_units": "nominal x/c,y/c", "ordering": "upper TE-LE-lower TE",
                               "source_urls": ["https://www.pdas.com/naca456thick4.html",
                                               "https://www.pdas.com/naca456mean2.html"],
                               "coefficient_last": -.1036 if closed else -.1015,
                               "reference_leading_edge_index": count // 2})
        return generated[key]

    profiles.extend([("Analítico NACA 0012", analytic("0012", 161)),
                     ("Analítico NACA 4412", analytic("4412", 161))])
    external: list[Airfoil] = []
    if fetch or reference_dir is not None:
        source = output / "originals" if fetch else reference_dir.resolve()
        records = published_sources(source, fetch=fetch)
        originals = output / "originals"
        originals.mkdir(exist_ok=True)
        for record in records:
            path = originals / record["file"]
            if not fetch:
                shutil.copyfile(source / record["file"], path)
            foil = Airfoil.from_dat(path, leading_edge_index=record["reference_leading_edge_index"])
            if len(foil.original_points) != record["point_count"]:
                raise AssertionError("published sample count changed")
            provenance.append({**record, "file": str(path.relative_to(output))})
            external.append(foil)
            profiles.append((f"UIUC {record['name']}", foil))
        _json(originals / "provenance.json", records)

    metrics = []
    fig, axes = plt.subplots((len(profiles) + 1) // 2, 2, figsize=(13, 3 * ((len(profiles) + 1)//2)),
                             squeeze=False, layout="constrained")
    for (label, foil), ax in zip(profiles, axes.flat):
        auto = Airfoil.from_points(foil.original_points)
        le_difference = float(np.linalg.norm(auto.leading_edge - foil.leading_edge) / foil.chord)
        metrics.append({"profile": label, "points": len(foil.original_points), "original_chord": foil.chord,
                        "area_normalized": foil.normalized_area, "perimeter_normalized": _perimeter(foil.normalized_points),
                        "explicit_le_index": foil.leading_edge_index, "auto_le_index": auto.leading_edge_index,
                        "auto_le_reference_difference_over_c": le_difference,
                        "auto_chord_relative_difference": auto.chord / foil.chord - 1,
                        "te_normalized_reference_error": float(np.linalg.norm(
                            (foil.normalized_points[0] + foil.normalized_points[-1]) / 2 - [1., 0.])),
                        "finite_te": not foil.is_te_closed})
        p, a = foil.normalized_points, auto.normalized_points
        ax.plot(p[:, 0], p[:, 1], label="LE explícita", linewidth=1.8)
        ax.plot(a[:, 0], a[:, 1], "--", label="LE automática", linewidth=1.)
        ax.scatter(0, 0, marker="x", color="black", label="Referência LE")
        ax.set_title(f"{label} — {len(p)} pontos\nDiferença LE/c = {le_difference:.3g}", fontsize=10)
        ax.set(xlabel="x/c próprio", ylabel="y/c próprio", aspect="equal")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.5), ncol=3)
        normalization.extend(_invariance(foil))
    for ax in list(axes.flat)[len(profiles):]:
        ax.set_visible(False)
    fig.suptitle("Referências geométricas — sem resultados aerodinâmicos")
    save_figure(fig, output / "normalized_profiles")
    plt.close(fig)
    _table(output, "reference_metrics", metrics)
    _table(output, "normalization_invariance", normalization)

    sensitivity, continuous = [], {}
    for code in ("0012", "4412"):
        for closed in (False, True):
            key = f"{code}_{'closed' if closed else 'finite'}"
            reference = continuous_metrics(code, closed)
            slot = continuous_slot_distance(code, closed)
            continuous[key] = {**reference, "slot": slot}
            for spacing in ("cosine", "uniform"):
                for count in (41, 81, 161, 321):
                    foil = analytic(code, count, closed, spacing)
                    main = analytic("0012", count, False, spacing)
                    geometry = MultiElementAirfoil(main, foil, FlapParameters(.3, (1.01, 0.)))
                    report = geometry.validate(min_clearance=.002)
                    auto = Airfoil.from_points(foil.original_points)
                    area, perimeter = foil.normalized_area, _perimeter(foil.normalized_points)
                    sensitivity.append({"code": code, "closed_te": closed, "spacing": spacing, "points": count,
                                        "area": area, "area_reference": reference["area"],
                                        "area_error": area - reference["area"],
                                        "perimeter": perimeter, "perimeter_reference": reference["perimeter"],
                                        "perimeter_error": perimeter - reference["perimeter"],
                                        "auto_le_index": auto.leading_edge_index,
                                        "explicit_le_index": foil.leading_edge_index,
                                        "auto_le_difference": float(np.linalg.norm(auto.leading_edge)),
                                        "te_midpoint_error": float(np.linalg.norm(foil.trailing_edge - [1., 0.])),
                                        "min_distance": report.min_distance, "continuous_distance": slot["distance"],
                                        "distance_error": report.min_distance - slot["distance"], "valid": report.valid})
    _table(output, "discretization", sensitivity)
    _json(output / "continuous_references.json", continuous)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout="constrained")
    for code in ("0012", "4412"):
        for spacing in ("cosine", "uniform"):
            rows = [r for r in sensitivity if r["code"] == code and not r["closed_te"] and r["spacing"] == spacing]
            for ax, field, title in zip(axes, ("area_error", "perimeter_error", "distance_error"),
                                        ("Erro absoluto de área / c²", "Erro absoluto de perímetro / c", "Erro absoluto de distância / c")):
                # Zero distance error for symmetric profiles is valid; floor only for display.
                ax.loglog([r["points"] for r in rows], [max(abs(r[field]), 1e-16) for r in rows],
                          "o-", label=f"{code} / {spacing}")
                ax.set(xlabel="Número de pontos", ylabel=title)
                ax.grid(which="both", alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle("Discretização vs referência contínua independente (piso gráfico: 1e-16)")
    save_figure(fig, output / "discretization_errors")
    plt.close(fig)

    cases = validation_cases(symmetric=analytic("0012", 161), cambered=analytic("4412", 81))
    if external:
        cases.extend([(f"published_{i+1}", MultiElementAirfoil(main, flap,
                      FlapParameters(.3, (1.1, -.04), angle)), .002, True, set())
                      for i, (main, flap, angle) in enumerate(zip(external, external[1:] + external[:1], (0.,25.,-20.), strict=True))])
    case_results, roundtrips = [], []
    for label, geometry, clearance, expected, expected_codes in cases:
        report = geometry.validate(min_clearance=clearance)
        codes = {issue.code for issue in report.issues}
        if report.valid != expected or not expected_codes.issubset(codes):
            raise AssertionError(f"Unexpected classification {label}: {report.to_dict()}")
        case_results.append({"id": label, "parameters": geometry.parameters.to_dict(),
                             "main_profile": geometry.main.name, "flap_profile": geometry.flap.name,
                             "validation": report.to_dict(), "expected_valid": expected})
        if report.valid:
            path = geometry.export(output / "exports" / label, configuration_id=label, min_clearance=clearance)
            _, main = read_dat(path / "main.dat")
            _, flap = read_dat(path / "flap.dat")
            np.testing.assert_array_equal(main, geometry.main_points)
            np.testing.assert_array_equal(flap, geometry.flap_points)
            recomputed = validate_elements(main, flap, min_clearance=clearance)
            meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
            if meta["validation"] != json.loads(json.dumps(recomputed.to_dict())):
                raise AssertionError("round-trip validation differs from metadata")
            if meta["parameters"] != json.loads(json.dumps(geometry.parameters.to_dict())):
                raise AssertionError("round-trip parameters differ from source geometry")
            roundtrips.append({"id": label, "main_max_point_error": float(np.max(abs(main - geometry.main_points))),
                               "flap_max_point_error": float(np.max(abs(flap - geometry.flap_points))),
                               "distance_error": abs(recomputed.min_distance - report.min_distance),
                               "finite_main_te": not geometry.main.is_te_closed,
                               "finite_flap_te": not geometry.flap.is_te_closed})
        print(f"{label:22} {'VALID' if report.valid else 'REJECTED':8} distance={report.min_distance}: {sorted(codes)}")
    _json(output / "configuration_results.json", case_results)
    _table(output, "roundtrips", roundtrips)
    for start in range(0, len(cases), 6):
        group = cases[start:start+6]
        fig = compare_configurations([row[1] for row in group], labels=[row[0] for row in group])
        # Each case has its own clearance: titles must use the SAME policy as its report.
        for ax, row in zip(fig.axes, group):
            label, geometry, clearance, _, _ = row
            state = "VÁLIDA" if geometry.validate(min_clearance=clearance).valid else "REJEITADA"
            parameters_title = ax.get_title().split("\n", 1)[1]
            ax.set_title(f"{label} — {state} (h/c={clearance:g})\n{parameters_title}", fontsize=10)
        save_figure(fig, output / f"configurations_{start//6+1}")
        plt.close(fig)
    _json(output / "provenance.json", provenance)
    summary = {"aeroopt_version": __version__, "python": platform.python_version(),
               "numpy": np.__version__, "shapely": shapely.__version__, "geos": shapely.geos_version_string,
               "matplotlib": matplotlib.__version__, "profile_count": len(profiles),
               "published_profile_count": len(external), "sensitivity_rows": len(sensitivity),
               "normalization_cases": len(normalization), "configurations": len(cases),
               "valid": sum(r["validation"]["valid"] for r in case_results), "roundtrips": len(roundtrips),
               "max_normalized_point_error": max(r["max_normalized_point_error"] for r in normalization),
               "max_relative_distance_error": max(r["max_relative_distance_error"] for r in normalization),
               "max_relative_area_error": max(r["relative_area_error"] for r in normalization),
               "license_policy": "published bytes and derivatives remain in ignored outputs; permission not inferred",
               "scope": "geometry only; no aerodynamic effectiveness established"}
    _json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    print(f"Results: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/validation_v011")
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--fetch-uiuc", action="store_true")
    sources.add_argument("--reference-dir", type=Path)
    args = parser.parse_args()
    run(args.output, args.reference_dir, args.fetch_uiuc)


if __name__ == "__main__":
    main()
