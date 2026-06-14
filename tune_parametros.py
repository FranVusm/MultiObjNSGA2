from __future__ import annotations

import argparse
import itertools
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


USED_FILE = Path("tuning_combinaciones_usadas.txt")
LOCK_FILE = Path("tuning_combinaciones_usadas.lock")
EXECUTABLE = Path("nsga2r_port.exe")


@dataclass(frozen=True)
class Combo:
    model: str
    popsize: int
    ngen: int
    pcross: float
    pmut: float

    @property
    def key(self) -> tuple[str, int, int, str, str]:
        return (
            self.model,
            self.popsize,
            self.ngen,
            format_float(self.pcross),
            format_float(self.pmut),
        )


def format_float(value: float) -> str:
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text if text else "0"


def model_name(model_path: Path) -> str:
    return model_path.stem


def acquire_lock() -> None:
    try:
        fd = os.open(str(LOCK_FILE), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(
            f"Ya existe {LOCK_FILE}. Parece que hay otro tuning corriendo. "
            "Si estas seguro de que no hay otro proceso, borra ese .lock."
        )
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(f"pid={os.getpid()} started={time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def release_lock() -> None:
    try:
        LOCK_FILE.unlink()
    except FileNotFoundError:
        pass


def ensure_used_file() -> None:
    if USED_FILE.exists():
        return
    USED_FILE.write_text(
        "# model,popsize,ngen,pcross,pmut,runs,status,feasible_count,result_dir,timestamp\n",
        encoding="utf-8",
    )


def read_used_keys() -> set[tuple[str, int, int, str, str]]:
    if not USED_FILE.exists():
        return set()

    used: set[tuple[str, int, int, str, str]] = set()
    for raw_line in USED_FILE.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",", 9)
        if len(parts) < 5:
            continue
        try:
            used.add((parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]))
        except ValueError:
            continue
    return used


def append_used(combo: Combo, runs: int, status: str, feasible_count: int | str, result_dir: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with USED_FILE.open("a", encoding="utf-8") as handle:
        handle.write(
            ",".join(
                [
                    combo.model,
                    str(combo.popsize),
                    str(combo.ngen),
                    format_float(combo.pcross),
                    format_float(combo.pmut),
                    str(runs),
                    status,
                    str(feasible_count),
                    result_dir,
                    timestamp,
                ]
            )
            + "\n"
        )


def candidate_combos(model: str) -> list[Combo]:
    priority_values = [
        (40, 100, 0.90, 0.20),
        (60, 100, 0.90, 0.20),
        (80, 100, 0.90, 0.20),
        (100, 100, 0.90, 0.20),
        (60, 150, 0.90, 0.20),
        (80, 150, 0.90, 0.20),
        (100, 150, 0.90, 0.20),
        (120, 150, 0.90, 0.20),
        (80, 200, 0.90, 0.20),
        (100, 200, 0.90, 0.20),
        (120, 200, 0.90, 0.20),
        (160, 200, 0.90, 0.20),
        (80, 250, 0.95, 0.20),
        (100, 250, 0.95, 0.20),
        (120, 250, 0.95, 0.20),
        (160, 250, 0.95, 0.20),
        (100, 300, 0.95, 0.15),
        (120, 300, 0.95, 0.15),
        (160, 300, 0.95, 0.15),
        (200, 300, 0.95, 0.15),
    ]

    combos: list[Combo] = [
        Combo(model, popsize, ngen, pcross, pmut)
        for popsize, ngen, pcross, pmut in priority_values
    ]

    grid = itertools.product(
        [40, 60, 80, 100, 120, 160, 200],
        [50, 100, 150, 200, 250, 300, 400],
        [0.80, 0.85, 0.90, 0.95],
        [0.10, 0.15, 0.20, 0.25, 0.30],
    )
    seen = {combo.key for combo in combos}
    for popsize, ngen, pcross, pmut in grid:
        combo = Combo(model, popsize, ngen, pcross, pmut)
        if combo.key not in seen:
            combos.append(combo)
            seen.add(combo.key)
    return combos


def parse_output_dir(output: str) -> str:
    matches = re.findall(r"Output directory:\s*(.+)", output)
    if matches:
        return matches[-1].strip()
    return ""


def parse_feasible_count(result_dir: str) -> int:
    if not result_dir:
        return 0
    report = Path(result_dir) / "aggregate_pareto_solutions.txt"
    if not report.exists():
        return 0
    for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("Feasible candidate solutions:"):
            return int(line.split(":", 1)[1].strip())
    return 0


def run_combo(combo: Combo, model_path: Path, seed: float, runs: int) -> tuple[int, str, int]:
    command = [
        str(EXECUTABLE),
        format_float(seed),
        str(model_path),
        str(combo.popsize),
        str(combo.ngen),
        format_float(combo.pcross),
        format_float(combo.pmut),
        str(runs),
    ]

    print("\nProbando combinacion:")
    print(
        f"  popsize={combo.popsize} ngen={combo.ngen} "
        f"pcross={format_float(combo.pcross)} pmut={format_float(combo.pmut)} runs={runs}"
    )
    print("  comando:", " ".join(command))

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    output_parts: list[str] = []
    assert process.stdout is not None
    for chunk in process.stdout:
        print(chunk, end="")
        output_parts.append(chunk)
    return_code = process.wait()

    output = "".join(output_parts)
    result_dir = parse_output_dir(output)
    feasible_count = parse_feasible_count(result_dir)
    return return_code, result_dir, feasible_count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prueba combinaciones de parametros hasta encontrar soluciones factibles."
    )
    parser.add_argument(
        "--model",
        default=r"..\..\data\Modelo_diverso.dat",
        help="Ruta del archivo .dat a probar.",
    )
    parser.add_argument("--seed", type=float, default=0.123, help="Seed base.")
    parser.add_argument("--runs", type=int, default=10, help="Seeds por combinacion.")
    parser.add_argument("--max-attempts", type=int, default=20, help="Intentos nuevos maximos.")
    parser.add_argument("--dry-run", action="store_true", help="Solo muestra las combinaciones nuevas.")
    args = parser.parse_args()

    if args.runs < 1:
        raise ValueError("--runs debe ser positivo")
    if args.max_attempts < 1:
        raise ValueError("--max-attempts debe ser positivo")
    if not EXECUTABLE.exists():
        raise FileNotFoundError(f"No existe {EXECUTABLE}. Compila nsga2r_port.exe primero.")

    model_path = Path(args.model)
    model = model_name(model_path)

    if args.dry_run:
        used = read_used_keys()
        available = [combo for combo in candidate_combos(model) if combo.key not in used]
        print(f"Combinaciones nuevas disponibles: {len(available)}")
        for combo in available[: args.max_attempts]:
            print(
                f"{combo.model}: popsize={combo.popsize}, ngen={combo.ngen}, "
                f"pcross={format_float(combo.pcross)}, pmut={format_float(combo.pmut)}"
            )
        return 0

    acquire_lock()
    try:
        ensure_used_file()
        used = read_used_keys()
        attempts = 0

        for combo in candidate_combos(model):
            if combo.key in used:
                continue
            if attempts >= args.max_attempts:
                break

            attempts += 1
            append_used(combo, args.runs, "RUNNING", "?", "")
            return_code, result_dir, feasible_count = run_combo(combo, model_path, args.seed, args.runs)

            status = "OK" if return_code == 0 else f"ERROR_{return_code}"
            if feasible_count > 0:
                status = "FOUND_FEASIBLE"
            append_used(combo, args.runs, status, feasible_count, result_dir)

            if feasible_count > 0:
                print("\nEncontrada combinacion con soluciones factibles:")
                print(
                    f"  popsize={combo.popsize} ngen={combo.ngen} "
                    f"pcross={format_float(combo.pcross)} pmut={format_float(combo.pmut)}"
                )
                print(f"  soluciones factibles: {feasible_count}")
                print(f"  resultados: {result_dir}")
                return 0

            print(f"\nSin factibles en esta combinacion. Intentos usados: {attempts}/{args.max_attempts}")

        print("\nNo se encontraron soluciones factibles en los intentos nuevos disponibles.")
        print(f"Intentos realizados en esta ejecucion: {attempts}")
        print(f"Historial: {USED_FILE}")
        return 1
    finally:
        release_lock()


if __name__ == "__main__":
    raise SystemExit(main())
