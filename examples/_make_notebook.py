#!/usr/bin/env python3
"""Write examples/Bell_examples.ipynb by running the teaching scripts.

From the repository root:

    python examples/_make_notebook.py

Each of examples/01_* ... examples/15_* contributes a markdown cell and a code
cell. The code cell includes that folder's run.jl. For 13_acceleration the
cell is examples/_notebook_13_smoke.jl, and 13_acceleration/run.jl is left as
the full ramp. Captured stdout and stderr are stored on the cells.

A successful run is cached under examples/_nb_cache/ so a later invocation
reruns only a script that changed or that failed. Set BELL_NOTEBOOK_FRESH=1
to ignore the cache.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
CACHE = EXAMPLES / "_nb_cache"
NOTEBOOK = EXAMPLES / "Bell_examples.ipynb"
SMOKE_NAME = "13_acceleration"
SMOKE_FILE = EXAMPLES / "_notebook_13_smoke.jl"

INTRO = {
    "01_reference_and_coordinates": (
        "Reference particle on `lattices/fodo.bell`, the longitudinal term "
        "`energy_term`, and the map between laboratory coordinates "
        "`(x, x', y, y', phi, delta)` and canonical coordinates "
        "`(x, px, y, py, l, delta)`. The conversion is also expanded on the "
        "coordinate jet, and the one-turn Taylor map is rewritten with "
        "`laboratory_map`."
    ),
    "02_lattice_io": (
        "Load `fodo.bell` and `magnetic_2.bell`, write the FODO line back to "
        "a structure file, insert a thin cavity, drop it with "
        "`magnetic_lattice`, and replace the reference particle."
    ),
    "03_element_hamiltonians": (
        "The s-Hamiltonian of a drift, a quadrupole, and a sector bend from "
        "the FODO file, plus a sextupole and a thin and a thick RF cavity "
        "built in the script. A thin cavity has no continuous Hamiltonian; "
        "its generator is `rf_potential`."
    ),
    "04_linear_optics": (
        "Linear optics of the 16-cell FODO line, then the same readings on "
        "`magnetic_2.bell`: element matrices, the one-turn matrix, periodic "
        "Twiss, dispersion, and the normal-form residual of that matrix."
    ),
    "05_element_maps": (
        "One element at a time. The Lie series of `H`, the Taylor map, the "
        "DA map, and `element_flow` / `exact_track_element` on one numeric "
        "point. A thin cavity is `rf_kick`."
    ),
    "06_lie_maps": (
        "Dragt-Finn factors of the FODO ring: `concatenate`, "
        "`product_factors`, and `track`. The degree-4 term of `concatenate` "
        "is checked against the Poisson bracket. `magnetic_2` then compares "
        "the Lie product with `ring_taylor_map` and a short tracked action."
    ),
    "07_map_readings": (
        "Readings of the one-turn map of `magnetic_2` with sextupoles at "
        "zero strength: second-order dispersion, momentum compaction, the "
        "slip factor, and chromaticity, each checked against tracking."
    ),
    "08_da_basics": (
        "Differential-algebra arithmetic in two variables, then "
        "`element_da_map` of a FODO quadrupole and `one_turn_map_da` of the "
        "16-cell line. `chromaticity_series` reads the order-3 map."
    ),
    "09_method_comparison": (
        "The same FODO elements, cell, and ring at orders 1-3. The Lie path, "
        "the DA map of `H`, the Taylor check, and numeric tracking. "
        "Dragt-Finn generators are taken from those maps."
    ),
    "10_sextupole_families": (
        "Place three sextupole families on `magnetic_2`, build the response "
        "in `k2`, solve two targets, and fit the same scalar residual with "
        "`NelderMead`, `RandomStart`, and `GridSearch`."
    ),
    "11_trackers_and_beam": (
        "Three particles in a `Beam` on the FODO ring. `FieldTracker` "
        "integrates the element s-Hamiltonian, `LieTracker` applies the "
        "stored factors, and `MapTracker` repeats the linear one-turn matrix."
    ),
    "12_rf_longitudinal": (
        "Fixed energy on `magnetic_2`. A thin cavity is `rf_kick`. A thick "
        "cavity is integrated from its s-Hamiltonian. The smooth model in "
        "`(phi, delta)` is the averaged comparison: `synchrotron_tune`, "
        "`long_H`, `separatrix`, and `bucket_height`."
    ),
    "13_acceleration": (
        "Thin-cavity acceleration. `13_acceleration/run.jl` is the full ramp "
        "from 10 MeV to 270 MeV, with an ensemble and the figures. The cell "
        "below is `examples/_notebook_13_smoke.jl`: the cavity check, 32 "
        "turns at fixed energy, and `accelerate` from 10 MeV to 10.3 MeV "
        "with two particles."
    ),
    "14_normal_form": (
        "One-turn Lie factor of `magnetic_2` with a thin cavity, then the "
        "homological normal form on the generators. Load the research module "
        "with `using Bell.NormalForm`."
    ),
    "15_fringe_fields": (
        "Fringe-field models on `magnetic_2`. The lattice file stays fixed. "
        "`FringeModel` selects `NoEdge`, `HardEdge`, or `SoftEdge` with an "
        "Enge, tanh, linear, or tabulated profile, and a splitting, Magnus, "
        "or Taylor integrator."
    ),
}


def julia_executable() -> str:
    return os.environ.get("JULIA", "julia")


def source_lines(text: str) -> list[str]:
    if not text.endswith("\n"):
        text += "\n"
    return text.splitlines(keepends=True)


def folders() -> list[Path]:
    found = sorted(
        p for p in EXAMPLES.iterdir()
        if p.is_dir() and len(p.name) > 3 and p.name[2] == "_" and p.name[:2].isdigit()
    )
    names = [p.name for p in found]
    expected = [f"{i:02d}_" for i in range(1, 16)]
    if len(found) != 15 or any(not n.startswith(prefix) for n, prefix in zip(names, expected)):
        raise SystemExit(f"expected folders 01-15, found {names}")
    return found


def script_for(folder: Path) -> Path:
    if folder.name == SMOKE_NAME:
        return SMOKE_FILE
    script = folder / "run.jl"
    if not script.is_file():
        raise SystemExit(f"missing {script}")
    return script


def cell_source(folder: Path, script: Path) -> str:
    if folder.name == SMOKE_NAME:
        return script.read_text(encoding="utf-8")
    rel = script.relative_to(EXAMPLES).as_posix()
    return (
        "include(joinpath(get(ENV, \"BELL_EXAMPLES_DIR\", pwd()), "
        f"\"{rel}\"))\n"
    )


def cache_path(name: str) -> Path:
    return CACHE / f"{name}.json"


def load_cache(name: str, script: Path) -> dict | None:
    if os.environ.get("BELL_NOTEBOOK_FRESH") == "1":
        return None
    path = cache_path(name)
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("mtime_ns") != script.stat().st_mtime_ns:
        return None
    if data.get("returncode") != 0:
        return None
    return data


def store_cache(name: str, script: Path, payload: dict) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    payload = dict(payload)
    payload["mtime_ns"] = script.stat().st_mtime_ns
    cache_path(name).write_text(json.dumps(payload), encoding="utf-8")


def run_script(name: str, script: Path) -> dict:
    cached = load_cache(name, script)
    if cached is not None:
        print(f"cache {name}  {cached['seconds']:.2f} s", flush=True)
        return cached
    env = os.environ.copy()
    env["BELL_EXAMPLES_DIR"] = str(EXAMPLES)
    env.setdefault("GKSwstype", "100")
    print(f"run   {name}  {script.name}", flush=True)
    started = time.perf_counter()
    proc = subprocess.run(
        [julia_executable(), f"--project={EXAMPLES}", str(script)],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    seconds = time.perf_counter() - started
    payload = {
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "seconds": seconds,
    }
    print(f"done  {name}  {seconds:.2f} s  exit {proc.returncode}", flush=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.stderr.write(proc.stdout[-4000:])
    else:
        store_cache(name, script, payload)
    return payload


def markdown_cell(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source_lines(text)}


def code_cell(source: str, stdout: str, stderr: str, count: int) -> dict:
    outputs = []
    if stdout:
        outputs.append({"output_type": "stream", "name": "stdout", "text": source_lines(stdout)})
    if stderr:
        outputs.append({"output_type": "stream", "name": "stderr", "text": source_lines(stderr)})
    return {
        "cell_type": "code",
        "execution_count": count,
        "metadata": {},
        "source": source_lines(source),
        "outputs": outputs,
    }


def main() -> int:
    chosen = folders()
    cells = [markdown_cell(
        "Bell examples, in teaching order.\n\n"
        "From the repository root, one script is "
        "`julia --project=examples examples/NN_name/run.jl`. "
        "This notebook was executed by `examples/_make_notebook.py`. "
        "A code cell includes that `run.jl`, except `13_acceleration`, "
        "whose cell is the short smoke in `examples/_notebook_13_smoke.jl`. "
        "The full ramp stays in `examples/13_acceleration/run.jl`.\n\n"
        "Re-running a cell needs the kernel working directory to be "
        "`examples/`, or `BELL_EXAMPLES_DIR` set to that directory.\n"
    )]
    count = 1
    failed = []
    timings = []
    for folder in chosen:
        script = script_for(folder)
        intro = INTRO.get(folder.name)
        if intro is None:
            raise SystemExit(f"no introduction for {folder.name}")
        rel = (folder / "run.jl").relative_to(ROOT).as_posix()
        cells.append(markdown_cell(f"## `{folder.name}`\n\n{intro}\n\n`{rel}`\n"))
        result = run_script(folder.name, script)
        timings.append((folder.name, result["seconds"], result["returncode"]))
        cells.append(code_cell(cell_source(folder, script), result["stdout"], result["stderr"], count))
        count += 1
        if result["returncode"] != 0:
            failed.append(folder.name)
    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Julia",
                "language": "julia",
                "name": "julia",
            },
            "language_info": {"name": "julia", "file_extension": ".jl"},
        },
        "cells": cells,
    }
    NOTEBOOK.write_text(json.dumps(notebook, indent=1), encoding="utf-8")
    print(f"wrote {NOTEBOOK}")
    for name, seconds, code in timings:
        print(f"  {name:32s} {seconds:8.2f} s  exit {code}")
    if failed:
        print("failed: " + ", ".join(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
