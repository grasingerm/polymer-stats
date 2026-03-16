#!/usr/bin/env python3
"""
Visualize the Hermans orientation parameter as a function of applied electric
field (E0) and transverse dipole susceptibility (K2) from MCMC simulation
output files.

The Hermans orientation parameter S is computed as:
    S = 0.5 * (3 * <cos²θ> / N_monomers - 1)

where <cos²θ> is the ensemble-averaged sum over all monomers (as reported
directly in the output file) and N_monomers = 100.

Usage:
    python plot_hermans.py <input_directory> [options]
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


# ---------------------------------------------------------------------------
# Parsing utilities
# ---------------------------------------------------------------------------

def parse_filename_params(filename: str) -> dict[str, float]:
    """Extract simulation parameters from the filename encoding.

    Convention: each parameter token is ``NAME-VALUE`` separated by
    underscores, where VALUE is an integer equal to 1000× the physical
    value.  The trailing ``run-N`` token is treated separately.
    """
    stem = os.path.splitext(os.path.basename(filename))[0]
    tokens = stem.split("_")
    params = {}
    for token in tokens:
        match = re.match(r"^([A-Za-z]\w*)-(\d+)$", token)
        if match:
            name, raw = match.group(1), match.group(2)
            if name == "run":
                params["run"] = int(raw)
            else:
                params[name] = int(raw) / 1000.0
    return params


def parse_output_file(filepath: str) -> dict[str, object]:
    """Read an MCMC .out file and return a dict of output quantities."""
    data = {}
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            # Try parsing as a Python literal (handles lists, floats, sci notation)
            try:
                data[key] = eval(value)  # safe here: controlled input files
            except Exception:
                data[key] = value
    return data


def hermans_parameter(cos2_theta_sum: float, n_monomers: int = 100) -> float:
    """Compute the Hermans orientation parameter from the monomer-summed
    ensemble average of cos²θ."""
    return 0.5 * (3.0 * cos2_theta_sum / n_monomers - 1.0)


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_data(input_dir: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scan *input_dir* for .out files and build the (E0, K2, S) grid.

    Returns
    -------
    E0_grid, K2_grid : 2-D arrays for meshgrid coordinates.
    S_grid           : 2-D array of Hermans parameter (averaged over runs).
    """
    # Accumulate S values keyed by (E0, K2)
    accumulator: dict[tuple[float, float], list[float]] = defaultdict(list)

    n_files = 0
    for fname in os.listdir(input_dir):
        if not fname.endswith(".out"):
            continue
        fpath = os.path.join(input_dir, fname)
        params = parse_filename_params(fname)
        if "E0" not in params or "K2" not in params:
            continue

        output = parse_output_file(fpath)
        cos2_key = "<cos2(θ)>"
        if cos2_key not in output:
            # Try a fallback for possible encoding differences
            cos2_key = next((k for k in output if "cos2" in k), None)
        if cos2_key is None:
            print(f"WARNING: no cos²θ entry in {fname}, skipping.",
                  file=sys.stderr)
            continue

        S = hermans_parameter(float(output[cos2_key]))
        accumulator[(params["E0"], params["K2"])].append(S)
        n_files += 1

    if n_files == 0:
        sys.exit(f"ERROR: no valid .out files found in '{input_dir}'.")

    # Build sorted unique axes
    E0_vals = sorted({key[0] for key in accumulator})
    K2_vals = sorted({key[1] for key in accumulator})

    E0_arr = np.array(E0_vals)
    K2_arr = np.array(K2_vals)
    E0_grid, K2_grid = np.meshgrid(E0_arr, K2_arr, indexing="ij")

    S_grid = np.full_like(E0_grid, np.nan)
    for i, e0 in enumerate(E0_vals):
        for j, k2 in enumerate(K2_vals):
            vals = accumulator.get((e0, k2))
            if vals:
                S_grid[i, j] = np.mean(vals)

    print(f"Loaded {n_files} files  →  "
          f"{len(E0_vals)} E0 × {len(K2_vals)} K2 grid  "
          f"({sum(len(v) for v in accumulator.values())} data points, "
          f"{len(accumulator)} unique grid cells)")

    return E0_grid, K2_grid, S_grid


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_hermans(E0: np.ndarray, K2: np.ndarray, S: np.ndarray,
                 output_path: str | None = None,
                 use_contour: bool = False,
                 log_x: bool = False,
                 log_y: bool = False,
                 xlim: tuple[float, float] | None = None,
                 ylim: tuple[float, float] | None = None) -> None:
    """Produce a filled-contour or pcolormesh plot of S(E0, K2)."""
    fig, ax = plt.subplots(figsize=(7, 5.5))

    if use_contour:
        n_levels = 30
        vmin = np.nanmin(S)
        vmax = np.nanmax(S)
        levels = np.linspace(vmin, vmax, n_levels)
        cf = ax.contourf(E0, K2, S, levels=levels, cmap="RdYlBu_r")
        cs = ax.contour(E0, K2, S, levels=levels[::3],
                        colors="k", linewidths=0.4, alpha=0.5)
        ax.clabel(cs, inline=True, fontsize=7, fmt="%.3f")
    else:
        cf = ax.pcolormesh(E0, K2, S, cmap="RdYlBu_r", shading="nearest")

    cbar = fig.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label(r"Hermans parameter $S$", fontsize=12)

    # --- Log scale --------------------------------------------------------
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    # --- Axis limits ------------------------------------------------------
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(r"$E_0$", fontsize=12)
    ax.set_ylabel(r"$\chi_{\perp}$", fontsize=12)
    #ax.set_title(r"Hermans orientation parameter $S(E_0,\,K_2)$",
    #             fontsize=13)

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Figure saved to {output_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot the Hermans orientation parameter from MCMC output "
                    "files over the (E0, K2) parameter space.")
    parser.add_argument("input_dir",
                        help="Directory containing .out simulation files.")
    parser.add_argument("-o", "--output",
                        help="Save figure to this path (e.g. hermans.png). "
                             "If omitted, an interactive window is shown.")
    parser.add_argument("--contour", action="store_true",
                        help="Use filled contours instead of a pixel heat map.")
    parser.add_argument("--logx", action="store_true",
                        help="Use logarithmic scale for the E0 axis.")
    parser.add_argument("--logy", action="store_true",
                        help="Use logarithmic scale for the K2 axis.")
    parser.add_argument("--loglog", action="store_true",
                        help="Use logarithmic scale for both axes.")
    parser.add_argument("--xlim", type=float, nargs=2, metavar=("XMIN", "XMAX"),
                        help="Set E0 axis limits (e.g. --xlim 0.01 1.0).")
    parser.add_argument("--ylim", type=float, nargs=2, metavar=("YMIN", "YMAX"),
                        help="Set K2 axis limits (e.g. --ylim 0.1 5.0).")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: '{args.input_dir}' is not a directory.")

    log_x = args.logx or args.loglog
    log_y = args.logy or args.loglog

    E0, K2, S = collect_data(args.input_dir)

    # Warn if log scale is requested but zeros are present in the data
    if log_x and np.any(E0 == 0):
        print("WARNING: E0 grid contains zero values; these will be "
              "invisible on a log scale. Consider using --xlim to exclude "
              "the origin.", file=sys.stderr)
    if log_y and np.any(K2 == 0):
        print("WARNING: K2 grid contains zero values; these will be "
              "invisible on a log scale. Consider using --ylim to exclude "
              "the origin.", file=sys.stderr)

    plot_hermans(E0, K2, S,
                 output_path=args.output,
                 use_contour=args.contour,
                 log_x=log_x,
                 log_y=log_y,
                 xlim=tuple(args.xlim) if args.xlim else None,
                 ylim=tuple(args.ylim) if args.ylim else None)


if __name__ == "__main__":
    main()
