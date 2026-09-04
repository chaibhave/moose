#!/usr/bin/env python3
"""Compare bicrystal GB profiles with Whipple's constant-source solution."""

import argparse
import csv
import math
from pathlib import Path

from scipy.integrate import quad


def whipple_interface_concentration(depth, time, d_bulk, d_gb, thickness):
    """Return Whipple's grain-side concentration at the GB interface."""
    diffusivity_ratio = d_gb / d_bulk
    eta = depth / math.sqrt(d_bulk * time)
    half_width = 0.5 * thickness
    beta = (diffusivity_ratio - 1.0) * half_width / math.sqrt(d_bulk * time)

    if eta == 0.0:
        return 1.0

    # sigma - 1 is localized on the scale beta when D_gb / D_bulk is
    # large. Integrating in q = (sigma - 1) / beta prevents an adaptive
    # quadrature over [1, Delta] from missing that narrow contribution.
    def transformed_integrand(q):
        sigma = 1.0 + beta * q
        if sigma >= diffusivity_ratio:
            return 0.0
        argument = 0.5 * math.sqrt(
            (diffusivity_ratio - 1.0) / (diffusivity_ratio - sigma)
        ) * ((sigma - 1.0) / beta)
        return beta * (
            math.exp(-(eta * eta) / (4.0 * sigma))
            * math.erfc(argument)
            / sigma**1.5
        )

    q_max = (diffusivity_ratio - 1.0) / beta
    integral, _ = quad(
        transformed_integrand,
        0.0,
        min(40.0, q_max),
        epsabs=2e-11,
        epsrel=2e-9,
        limit=300,
    )
    return math.erfc(0.5 * eta) + eta * integral / (2.0 * math.sqrt(math.pi))


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="bicrystal_whipple")
    parser.add_argument("--d-bulk", type=float, default=1e-19)
    parser.add_argument("--d-gb", type=float, default=1e-13)
    parser.add_argument("--thickness", type=float, default=0.5e-9)
    # The checked-in 1 um TRI3 mesh gives 5.9e-4 or less at the four sample
    # times, leaving a factor-of-three margin for supported solver platforms.
    parser.add_argument("--relative-l2-tolerance", type=float, default=0.002)
    args = parser.parse_args()

    prefix = Path(args.prefix)
    time_rows = read_csv(prefix.with_name(prefix.name + "_gb_profile_time.csv"))
    results = []

    for time_row in time_rows:
        time = float(time_row["time"])
        step = int(time_row["timestep"])
        if time == 0.0:
            continue

        profile_path = prefix.with_name(prefix.name + f"_gb_profile_{step:04d}.csv")
        profile_rows = read_csv(profile_path)
        squared_error = 0.0
        squared_reference = 0.0
        maximum_error = 0.0
        comparison_rows = []

        for row in profile_rows:
            depth = float(row["x"])
            fem_value = float(row["c"])
            reference = whipple_interface_concentration(
                depth, time, args.d_bulk, args.d_gb, args.thickness
            )
            error = fem_value - reference
            squared_error += error * error
            squared_reference += reference * reference
            maximum_error = max(maximum_error, abs(error))
            comparison_rows.append((depth, fem_value, reference, error))

        relative_l2 = math.sqrt(squared_error / squared_reference)
        results.append((time, relative_l2, maximum_error))
        comparison_path = prefix.with_name(prefix.name + f"_whipple_{step:04d}.csv")
        with comparison_path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "fem", "whipple", "error"))
            writer.writerows(comparison_rows)

    summary_path = prefix.with_name(prefix.name + "_whipple_comparison.csv")
    with summary_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time", "relative_l2_error", "maximum_absolute_error"))
        writer.writerows(results)

    worst_error = max(row[1] for row in results)
    print(f"Whipple comparison: worst relative L2 error = {worst_error:.6g}")
    if worst_error > args.relative_l2_tolerance:
        raise SystemExit(
            f"relative L2 error exceeds tolerance {args.relative_l2_tolerance:g}"
        )


if __name__ == "__main__":
    main()
