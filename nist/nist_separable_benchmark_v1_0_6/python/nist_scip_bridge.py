"""PySCIPOpt bridge for bounded separable NIST StRD regression problems."""
from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Any

import numpy as np
from pyscipopt import Model, exp, quicksum

BRIDGE_API_VERSION = "1.0.0"
BRIDGE_MODULE_NAME = "nist_scip_bridge"


def _json(value: Any) -> str:
    return json.dumps(value, allow_nan=False, separators=(",", ":"))


def bridge_api_version() -> str:
    return BRIDGE_API_VERSION


def bridge_info_json() -> str:
    return _json({
        "api_version": BRIDGE_API_VERSION,
        "module_name": BRIDGE_MODULE_NAME,
        "module_file": str(Path(__file__).resolve()),
        "python_executable": sys.executable,
    })


def _package_version(name: str) -> str | None:
    try:
        return package_version(name)
    except PackageNotFoundError:
        return None


def _safe(obj: Any, name: str, converter=None):
    try:
        value = getattr(obj, name)()
        return converter(value) if converter else value
    except Exception:
        return None


def versions_json() -> str:
    model = Model("version_check")
    try:
        scip = f"{model.getMajorVersion()}.{model.getMinorVersion()}.{model.getTechVersion()}"
    finally:
        try:
            model.free()
        except Exception:
            pass
    return _json({
        "bridge_api_version": BRIDGE_API_VERSION,
        "bridge_module_name": BRIDGE_MODULE_NAME,
        "bridge_file": str(Path(__file__).resolve()),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "machine": platform.machine(),
        "platform": platform.platform(),
        "numpy": _package_version("numpy"),
        "pyscipopt": _package_version("pyscipopt"),
        "scip": scip,
        "thread_environment": {name: os.environ.get(name) for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
    })


def smoke_test_json() -> str:
    model = Model("nist_bridge_smoke_test")
    try:
        model.setParam("display/verblevel", 0)
        x = model.addVar(name="x", lb=0.0, ub=2.0)
        y = model.addVar(name="y", lb=0.0, ub=2.0)
        model.addCons(x + y >= 1.0)
        model.setObjective(x + y, sense="minimize")
        model.optimize()
        sol = model.getBestSol()
        return _json({
            "status": str(model.getStatus()),
            "has_solution": sol is not None,
            "objective": None if sol is None else float(model.getSolObjVal(sol)),
            "dual_bound": float(model.getDualbound()),
            "nodes": int(model.getNNodes()),
            "solve_time": float(model.getSolvingTime()),
        })
    finally:
        try:
            model.free()
        except Exception:
            pass


def _options(text: Any) -> dict[str, Any]:
    if text is None or not str(text):
        return {}
    value = json.loads(str(text))
    if not isinstance(value, dict):
        raise ValueError("options_json must encode an object")
    return value


def _vector(value: Any, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} is empty or non-finite")
    return result


def _matrix(value: Any, rows: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim == 1:
        if result.size % rows:
            raise ValueError(f"{name} has incompatible dimensions")
        result = np.reshape(result, (rows, result.size // rows), order="F")
    if result.ndim != 2 or result.shape[0] != rows or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} has incompatible dimensions")
    return result


def _evaluate(problem_id: str, X: np.ndarray, response: np.ndarray,
              linear: np.ndarray, z: np.ndarray, residual_scale: float) -> tuple[float, np.ndarray, np.ndarray]:
    x = X[:, 0]
    pid = problem_id.lower()
    if pid == "boxbod":
        prediction = linear[0] * (1.0 - np.exp(-z[0] * x))
    elif pid == "mgh17":
        prediction = linear[0] + linear[1] * np.exp(-z[0] * x) + linear[2] * np.exp(-z[1] * x)
    elif pid == "lanczos2":
        prediction = sum(linear[j] * np.exp(-z[j] * x) for j in range(3))
    elif pid == "kirby2":
        prediction = (linear[0] + linear[1] * x + linear[2] * x * x) / (1.0 + z[0] * x + z[1] * x * x)
    elif pid == "hahn1":
        prediction = (linear[0] + linear[1] * x + linear[2] * x**2 + linear[3] * x**3) / (1.0 + z[0] * x + z[1] * x**2 + z[2] * x**3)
    elif pid == "gauss1":
        prediction = (linear[0] * np.exp(-z[0] * x)
                      + linear[1] * np.exp(-((x - z[1]) / z[2]) ** 2)
                      + linear[2] * np.exp(-((x - z[3]) / z[4]) ** 2))
    elif pid == "nelson":
        x1, x2 = X[:, 0], X[:, 1]
        prediction = linear[0] - linear[1] * x1 * np.exp(-z[0] * x2)
    elif pid == "mgh10":
        prediction = linear[0] * np.exp(z[0] / (x + z[1]))
    else:
        raise ValueError(f"Unsupported NIST problem {problem_id}")
    prediction = prediction / residual_scale
    residual = prediction - response
    return float(residual @ residual), prediction, residual


def _prediction_expression(problem_id: str, row: np.ndarray, c: list[Any], z: list[Any], residual_scale: float):
    pid = problem_id.lower()
    x = float(row[0])
    if pid == "boxbod":
        value = c[0] * (1.0 - exp(-z[0] * x))
    elif pid == "mgh17":
        value = c[0] + c[1] * exp(-z[0] * x) + c[2] * exp(-z[1] * x)
    elif pid == "lanczos2":
        value = quicksum(c[j] * exp(-z[j] * x) for j in range(3))
    elif pid == "kirby2":
        value = (c[0] + c[1] * x + c[2] * x * x) / (1.0 + z[0] * x + z[1] * x * x)
    elif pid == "hahn1":
        value = (c[0] + c[1] * x + c[2] * x**2 + c[3] * x**3) / (1.0 + z[0] * x + z[1] * x**2 + z[2] * x**3)
    elif pid == "gauss1":
        value = (c[0] * exp(-z[0] * x)
                 + c[1] * exp(-((x - z[1]) / z[2]) ** 2)
                 + c[2] * exp(-((x - z[3]) / z[4]) ** 2))
    elif pid == "nelson":
        x1, x2 = float(row[0]), float(row[1])
        value = c[0] - c[1] * x1 * exp(-z[0] * x2)
    elif pid == "mgh10":
        value = c[0] * exp(z[0] / (x + z[1]))
    else:
        raise ValueError(f"Unsupported NIST problem {problem_id}")
    return value / residual_scale


def _configure(model: Model, options: dict[str, Any]) -> dict[str, Any]:
    time_limit = float(options.get("time_limit", 600.0))
    eps_r = float(options.get("paper_epsilon_relative", 0.01))
    eps_a = float(options.get("paper_epsilon_absolute", 1e-8))
    threads = int(options.get("threads", 1))
    verbosity = int(options.get("verbosity", 2))
    feastol = float(options.get("feasibility_tolerance", 1e-6))
    disable_presolve = bool(options.get("disable_presolve", True))
    native_gap = eps_r / (1.0 - eps_r) if eps_r > 0 else 0.0
    model.setParam("limits/time", time_limit)
    model.setParam("limits/gap", native_gap)
    model.setParam("limits/absgap", eps_a)
    model.setParam("display/verblevel", verbosity)
    model.setParam("parallel/minnthreads", threads)
    model.setParam("parallel/maxnthreads", threads)
    model.setParam("parallel/mode", 1)
    model.setParam("timing/clocktype", 2)
    model.setParam("numerics/feastol", feastol)
    model.setParam("randomization/randomseedshift", int(options.get("random_seed_shift", 0)))
    model.setParam("randomization/permutationseed", int(options.get("permutation_seed", 0)))
    model.setParam("randomization/lpseed", int(options.get("lp_seed", 0)))
    if disable_presolve:
        model.setParam("presolving/maxrounds", 0)
        model.setParam("presolving/maxrestarts", 0)
    return {
        "time_limit": time_limit,
        "paper_epsilon_relative": eps_r,
        "paper_epsilon_absolute": eps_a,
        "native_relative_gap_limit": native_gap,
        "native_absolute_gap_limit": eps_a,
        "threads": threads,
        "verbosity": verbosity,
        "feasibility_tolerance": feastol,
        "disable_presolve": disable_presolve,
    }


def solve_nist_json(problem_id_value: Any, predictors_value: Any, response_value: Any,
                    residual_scale_value: Any,
                    linear_lower_value: Any, linear_upper_value: Any,
                    z_lower_value: Any, z_upper_value: Any,
                    initial_linear_value: Any, initial_z_value: Any,
                    options_json: Any = None) -> str:
    total_start = time.perf_counter()
    problem_id = str(problem_id_value)
    response = _vector(response_value, "response")
    residual_scale = float(residual_scale_value)
    if not math.isfinite(residual_scale) or residual_scale <= 0:
        raise ValueError("residual_scale must be positive and finite")
    X = _matrix(predictors_value, response.size, "predictors")
    linear_lower = _vector(linear_lower_value, "linear_lower")
    linear_upper = _vector(linear_upper_value, "linear_upper")
    z_lower = _vector(z_lower_value, "z_lower")
    z_upper = _vector(z_upper_value, "z_upper")
    initial_linear = _vector(initial_linear_value, "initial_linear")
    initial_z = _vector(initial_z_value, "initial_z")
    if linear_lower.size != linear_upper.size or np.any(linear_lower >= linear_upper):
        raise ValueError("Invalid linear bounds")
    if z_lower.size != z_upper.size or np.any(z_lower >= z_upper):
        raise ValueError("Invalid nonlinear bounds")
    if initial_linear.size != linear_lower.size or initial_z.size != z_lower.size:
        raise ValueError("Invalid initial point dimensions")
    initial_linear = np.minimum(np.maximum(initial_linear, linear_lower), linear_upper)
    initial_z = np.minimum(np.maximum(initial_z, z_lower), z_upper)
    external_upper, initial_prediction, initial_residual = _evaluate(
        problem_id, X, response, initial_linear, initial_z, residual_scale)
    options = _options(options_json)

    model = Model(str(options.get("model_name", f"nist_{problem_id}")))
    settings = _configure(model, options)
    if options.get("log_file"):
        path = Path(str(options["log_file"])).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        model.setLogfile(str(path))

    c_vars = [model.addVar(name=f"c_{j+1}", lb=float(linear_lower[j]), ub=float(linear_upper[j]))
              for j in range(linear_lower.size)]
    z_vars = [model.addVar(name=f"z_{j+1}", lb=float(z_lower[j]), ub=float(z_upper[j]))
              for j in range(z_lower.size)]
    t_vars = []
    residual_exprs = []
    for i in range(response.size):
        prediction = _prediction_expression(problem_id, X[i, :], c_vars, z_vars, residual_scale)
        residual = prediction - float(response[i])
        residual_exprs.append(residual)
        t = model.addVar(name=f"sq_residual_{i+1}", lb=0.0)
        model.addCons(t >= residual * residual, name=f"rss_epigraph_{i+1}")
        t_vars.append(t)
    objective_expr = quicksum(t_vars)
    model.setObjective(objective_expr, sense="minimize")
    cutoff_margin = max(settings["paper_epsilon_absolute"],
                        10.0 * settings["feasibility_tolerance"] * max(1.0, external_upper))
    model.addCons(objective_expr <= external_upper + cutoff_margin, name="external_incumbent_cutoff")

    candidate = model.createOrigSol()
    for var, value in zip(c_vars, initial_linear):
        model.setSolVal(candidate, var, float(value))
    for var, value in zip(z_vars, initial_z):
        model.setSolVal(candidate, var, float(value))
    for var, value in zip(t_vars, initial_residual * initial_residual):
        model.setSolVal(candidate, var, float(value))
    candidate_stored = bool(model.addSol(candidate, free=True))

    if options.get("write_problem"):
        path = Path(str(options["write_problem"])).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        model.writeProblem(str(path), trans=False, genericnames=False, verbose=False)

    build_time = time.perf_counter() - total_start
    model.optimize()
    raw_status = str(model.getStatus())
    solution = model.getBestSol()

    selected = {
        "source": "external_initial", "objective": external_upper,
        "linear": initial_linear, "z": initial_z,
        "prediction": initial_prediction, "residual": initial_residual,
        "epigraph_objective": external_upper,
    }
    scip_candidate = None
    if solution is not None:
        c_raw = np.asarray([model.getSolVal(solution, v) for v in c_vars], dtype=float)
        z_raw = np.asarray([model.getSolVal(solution, v) for v in z_vars], dtype=float)
        c_projected = np.minimum(np.maximum(c_raw, linear_lower), linear_upper)
        z_projected = np.minimum(np.maximum(z_raw, z_lower), z_upper)
        objective, prediction, residual = _evaluate(
            problem_id, X, response, c_projected, z_projected, residual_scale)
        scip_candidate = {
            "source": "scip_recomputed_projected", "objective": objective,
            "linear": c_projected, "z": z_projected,
            "prediction": prediction, "residual": residual,
            "epigraph_objective": float(model.getSolObjVal(solution)),
            "linear_projection_inf_norm": float(np.max(np.abs(c_projected-c_raw))),
            "z_projection_inf_norm": float(np.max(np.abs(z_projected-z_raw))),
        }
        if objective < selected["objective"]:
            selected = scip_candidate

    raw_dual = float(model.getDualbound()) if math.isfinite(float(model.getDualbound())) else None
    upper = float(selected["objective"])
    dual = max(0.0, raw_dual if raw_dual is not None else 0.0)
    consistency_tol = max(1e-8, 10.0 * settings["feasibility_tolerance"] * max(1.0, abs(upper)))
    consistent = dual <= upper + consistency_tol
    margin = dual - ((1.0-settings["paper_epsilon_relative"])*upper-settings["paper_epsilon_absolute"])
    certificate = bool(consistent and margin >= -consistency_tol)
    if not consistent:
        status = "invalid_dual_bound"
    elif raw_status == "optimal":
        status = "optimal"
    elif raw_status in ("gaplimit", "absgaplimit") and certificate:
        status = "paper_gap_native"
    else:
        status = raw_status

    gap_value = _safe(model, "getGap", float)
    if gap_value is not None and not math.isfinite(gap_value):
        gap_value = None
    result = {
        "bridge_api_version": BRIDGE_API_VERSION,
        "bridge_module_file": str(Path(__file__).resolve()),
        "status": status,
        "scip_status": raw_status,
        "has_solution": True,
        "solution_source": selected["source"],
        "objective": upper,
        "primal_bound": upper,
        "dual_bound": dual,
        "raw_scip_dual_bound": raw_dual,
        "bound_consistent": consistent,
        "bound_consistency_tolerance": consistency_tol,
        "paper_certificate_satisfied": certificate,
        "paper_certificate_margin": margin,
        "scip_relative_gap": gap_value,
        "linear": selected["linear"].tolist(),
        "z": selected["z"].tolist(),
        "prediction": selected["prediction"].tolist(),
        "residual": selected["residual"].tolist(),
        "external_initial_objective": external_upper,
        "candidate_stored": candidate_stored,
        "scip_epigraph_objective": None if scip_candidate is None else scip_candidate["epigraph_objective"],
        "scip_objective_recomputed": None if scip_candidate is None else scip_candidate["objective"],
        "nodes": _safe(model, "getNNodes", int),
        "total_nodes": _safe(model, "getNTotalNodes", int),
        "solve_time": _safe(model, "getSolvingTime", float),
        "build_time": build_time,
        "settings": settings,
        "versions": json.loads(versions_json()),
        "model": {
            "problem_id": problem_id,
            "n_linear": int(linear_lower.size),
            "n_nonlinear": int(z_lower.size),
            "n_residual": int(response.size),
            "residual_scale": residual_scale,
        },
    }
    result["python_total_time"] = time.perf_counter() - total_start
    try:
        return _json(result)
    finally:
        try:
            model.free()
        except Exception:
            pass
