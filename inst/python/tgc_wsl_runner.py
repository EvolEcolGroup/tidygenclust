"""File-based bridge for running tidygenclust Python backends in WSL."""

import argparse
import json
import os
from pathlib import Path

def _write_matrix(output_dir, name, value):
    """Write a matrix as row-major little-endian float64 data."""
    import numpy as np

    matrix = np.asarray(value, dtype="<f8", order="C")
    if matrix.ndim != 2:
        raise ValueError(f"{name} result is not a two-dimensional matrix")
    filename = f"{name}.bin"
    matrix.tofile(output_dir / filename)
    return {"file": filename, "shape": list(matrix.shape)}


def main():
    # R writes one job file so paths and analysis options need no shell parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parsed = parser.parse_args()

    with parsed.config.open("r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    arguments = config["args"]
    namespace = argparse.Namespace(**arguments)

    # Set thread limits before either backend imports numerical libraries.
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = str(arguments["threads"])

    # Reuse the original wrappers rather than duplicate either analysis here.
    algorithm = config["algorithm"]
    if algorithm == "fastmixture":
        from py_rfastmixture import fastmixture_run

        result = fastmixture_run(namespace)
    elif algorithm == "adamixture":
        from py_radamixture import adamixture_run

        result = adamixture_run(args=namespace)
    else:
        raise ValueError(f"Unknown tidygenclust algorithm: {algorithm}")

    has_cv = arguments.get("cv") is not None
    # Match the return conventions used by the existing R/reticulate wrappers.
    p_matrix = cv_result = None
    if arguments["no_freqs"]:
        if has_cv:
            q_matrix, cv_result = result
        else:
            q_matrix = result
    elif has_cv:
        q_matrix, p_matrix, cv_result = result
    else:
        q_matrix, p_matrix = result

    # Keep large matrices binary; JSON records their filenames, shapes and CV.
    # R reads this manifest only after the analysis process exits successfully.
    manifest = {
        "Q": _write_matrix(output_dir, "Q", q_matrix),
        "P": None,
        "cv": None,
    }
    if p_matrix is not None:
        manifest["P"] = _write_matrix(output_dir, "P", p_matrix)
    if cv_result is not None:
        average = cv_result["avg"] if isinstance(cv_result, dict) else cv_result.avg
        manifest["cv"] = float(average)

    with (output_dir / "result.json").open("w", encoding="utf-8") as result_file:
        json.dump(manifest, result_file, allow_nan=False)


if __name__ == "__main__":
    main()
