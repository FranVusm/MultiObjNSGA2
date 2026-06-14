from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt


@dataclass
class ParetoPoint:
    index: int
    f1: float
    f2: float
    source_run: int
    seed: float
    feasible: bool


def latest_results_dir(results_root: Path) -> Path:
    if not results_root.exists():
        raise FileNotFoundError(f"No existe la carpeta de resultados: {results_root}")
    candidates = [p for p in results_root.iterdir() if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No hay subcarpetas de resultados en: {results_root}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def resolve_points_file(path_arg: str | None) -> Path:
    if path_arg is None:
        result_dir = latest_results_dir(Path("resultados"))
        return result_dir / "aggregate_pareto_points.csv"

    path = Path(path_arg)
    if path.is_dir():
        return path / "aggregate_pareto_points.csv"
    return path


def read_pareto_points(points_file: Path) -> list[ParetoPoint]:
    if not points_file.exists():
        raise FileNotFoundError(f"No existe el archivo: {points_file}")

    points: list[ParetoPoint] = []
    with points_file.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if not header.startswith("index,F1,F2"):
            raise ValueError(f"Archivo inesperado: {points_file}")

        for line_number, line in enumerate(handle, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split(",", 6)
            if len(fields) < 6:
                raise ValueError(f"Linea incompleta {line_number}: {line}")

            points.append(
                ParetoPoint(
                    index=int(fields[0]),
                    f1=float(fields[1]),
                    f2=float(fields[2]),
                    source_run=int(fields[3]),
                    seed=float(fields[4]),
                    feasible=fields[5].strip().lower() == "true",
                )
            )

    return sorted(points, key=lambda p: (p.f1, p.f2))


def plot_pareto(points: list[ParetoPoint], output_file: Path, title: str) -> None:
    if not points:
        raise ValueError("No hay puntos de Pareto para graficar.")

    feasible = [p for p in points if p.feasible]
    infeasible = [p for p in points if not p.feasible]

    fig, ax = plt.subplots(figsize=(9, 6), dpi=140)

    if len(points) > 1:
        ax.plot(
            [p.f1 for p in points],
            [p.f2 for p in points],
            color="#2f6f9f",
            linewidth=1.4,
            alpha=0.75,
            label="Frente agregado",
        )

    if feasible:
        ax.scatter(
            [p.f1 for p in feasible],
            [p.f2 for p in feasible],
            s=60,
            color="#1b9e77",
            edgecolor="black",
            linewidth=0.5,
            label="Factible",
            zorder=3,
        )

    if infeasible:
        ax.scatter(
            [p.f1 for p in infeasible],
            [p.f2 for p in infeasible],
            s=60,
            color="#d95f02",
            edgecolor="black",
            linewidth=0.5,
            label="No factible",
            zorder=3,
        )

    for point in points:
        ax.annotate(
            str(point.index),
            (point.f1, point.f2),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
        )

    ax.set_title(title)
    ax.set_xlabel("F1")
    ax.set_ylabel("F2")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.45)
    ax.legend()
    ax.margins(0.08)
    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Grafica el frente Pareto desde aggregate_pareto_points.csv."
    )
    parser.add_argument(
        "path",
        nargs="?",
        help="Carpeta de resultados o archivo aggregate_pareto_points.csv. Si se omite, usa la ultima carpeta en ./resultados.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Ruta de salida del grafico. Por defecto: pareto_front.png en la carpeta del CSV.",
    )
    parser.add_argument(
        "--title",
        default="Frente Pareto agregado",
        help="Titulo del grafico.",
    )
    args = parser.parse_args()

    points_file = resolve_points_file(args.path)
    points = read_pareto_points(points_file)
    output_file = Path(args.output) if args.output else points_file.with_name("pareto_front.png")

    plot_pareto(points, output_file, args.title)
    print(f"Grafico guardado en: {output_file}")
    print(f"Puntos graficados: {len(points)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
