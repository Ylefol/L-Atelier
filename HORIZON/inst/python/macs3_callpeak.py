"""
macs3_callpeak.py
-----------------
Python helper for HORIZON_call_peaks().

Arguments are injected from R via reticulate::py_set_attr() into the
module-level `args` dict before this script is run via py_run_file().

Expected keys in `args`:
    treatment   : str  — path to treatment BED
    control     : str or None — path to control BED (CHIP), else None
    output_dir  : str  — output directory
    sample_name : str  — prefix for output files
    genome_size : str or float — effective genome size ('hs', 'mm', or numeric)
    q_value     : float — FDR threshold (e.g. 0.05)
    extra_flags : list[str] — assay-type-specific flags
    broad       : bool — call broad peaks

Returns (set on module):
    result_peak : str — path to the narrowPeak or broadPeak file
"""

import subprocess
import sys
import os


def _run(args_dict):
    cmd = [sys.executable, "-m", "MACS3", "callpeak"]

    cmd += ["-t", args_dict["treatment"]]
    if args_dict.get("control"):
        cmd += ["-c", args_dict["control"]]

    cmd += [
        "--outdir",  args_dict["output_dir"],
        "-n",        args_dict["sample_name"],
        "-g",        str(args_dict["genome_size"]),
        "-q",        str(args_dict["q_value"]),
        "--keep-dup", "all",
        "-f",        "BED",
    ]

    cmd += list(args_dict.get("extra_flags", []))

    if args_dict.get("broad"):
        cmd += ["--broad", "--broad-cutoff", str(args_dict["q_value"])]

    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"MACS3 callpeak failed (exit code {result.returncode})"
        )

    ext = "broadPeak" if args_dict.get("broad") else "narrowPeak"
    peak_file = os.path.join(
        args_dict["output_dir"],
        f"{args_dict['sample_name']}_peaks.{ext}"
    )
    return peak_file


# `args` is injected by R before py_run_file; fall back to {} if running
# directly (e.g. for testing)
if "args" not in dir():
    args = {}

result_peak = _run(args)
