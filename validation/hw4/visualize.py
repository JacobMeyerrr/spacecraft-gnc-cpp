#!/usr/bin/env python3

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


# ================================================================
# ASTE 585 - HOMEWORK 4 visualization
#
# The C++ program does the numerical simulation.
# This Python script only reads the CSV output and makes the plots.
#
# Usage:
#     python3 validation/hw4/visualize.py design_a
#     python3 validation/hw4/visualize.py design_b
# ================================================================


ROOT = Path(__file__).resolve().parents[2]
HW4 = ROOT / "data" / "validation" / "hw4"


def load_csv(name):
    path = HW4 / name

    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run ./build/validate_hw4 first."
        )

    return np.genfromtxt(
        path,
        delimiter=",",
        names=True,
    )


def plot_step_response(prefix):
    x = load_csv(f"hw4_{prefix}_roll_step.csv")
    y = load_csv(f"hw4_{prefix}_yaw_step.csv")
    p = load_csv(f"hw4_{prefix}_pitch_step.csv")

    # X-axis step: roll + coupled yaw.
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
    )

    axes[0].plot(
        x["time_s"] / 3600.0,
        x["phi_open_deg"],
        label="Open loop",
    )

    axes[0].plot(
        x["time_s"] / 3600.0,
        x["phi_closed_deg"],
        "--",
        label="Closed loop",
    )

    axes[0].axhline(0.06, linestyle=":")
    axes[0].axhline(-0.06, linestyle=":")
    axes[0].set_ylabel("Roll phi [deg]")
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(
        x["time_s"] / 3600.0,
        x["psi_open_deg"],
        label="Open loop",
    )

    axes[1].plot(
        x["time_s"] / 3600.0,
        x["psi_closed_deg"],
        "--",
        label="Closed loop",
    )

    axes[1].axhline(0.5, linestyle=":")
    axes[1].axhline(-0.5, linestyle=":")
    axes[1].set_xlabel("Time [hr]")
    axes[1].set_ylabel("Yaw psi [deg]")
    axes[1].grid(True)
    axes[1].legend()

    fig.suptitle(
        f"HW4 {prefix}: X-axis Step Disturbance"
    )

    fig.tight_layout()

    # Z-axis step: yaw + coupled roll.
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
    )

    axes[0].plot(
        y["time_s"] / 3600.0,
        y["phi_open_deg"],
        label="Open loop",
    )

    axes[0].plot(
        y["time_s"] / 3600.0,
        y["phi_closed_deg"],
        "--",
        label="Closed loop",
    )

    axes[0].axhline(0.06, linestyle=":")
    axes[0].axhline(-0.06, linestyle=":")
    axes[0].set_ylabel("Roll phi [deg]")
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(
        y["time_s"] / 3600.0,
        y["psi_open_deg"],
        label="Open loop",
    )

    axes[1].plot(
        y["time_s"] / 3600.0,
        y["psi_closed_deg"],
        "--",
        label="Closed loop",
    )

    axes[1].axhline(0.5, linestyle=":")
    axes[1].axhline(-0.5, linestyle=":")
    axes[1].set_xlabel("Time [hr]")
    axes[1].set_ylabel("Yaw psi [deg]")
    axes[1].grid(True)
    axes[1].legend()

    fig.suptitle(
        f"HW4 {prefix}: Z-axis Step Disturbance"
    )

    fig.tight_layout()

    # Pitch step.
    fig, axis = plt.subplots(
        figsize=(9, 5),
    )

    axis.plot(
        p["time_s"],
        p["theta_open_deg"],
        label="Open loop",
    )

    axis.plot(
        p["time_s"],
        p["theta_closed_deg"],
        "--",
        label="Closed loop",
    )

    axis.axhline(0.06, linestyle=":")
    axis.axhline(-0.06, linestyle=":")
    axis.set_xlabel("Time [s]")
    axis.set_ylabel("Pitch theta [deg]")
    axis.set_title(
        f"HW4 {prefix}: Y-axis Step Disturbance"
    )
    axis.grid(True)
    axis.legend()

    fig.tight_layout()


def plot_periodic(prefix):
    data = load_csv(
        f"hw4_{prefix}_periodic.csv"
    )

    t_hr = data["time_s"] / 3600.0

    # Open vs closed attitude.
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(9, 8),
        sharex=True,
    )

    pairs = [
        ("phi_open_deg", "phi_closed_deg", "Roll phi [deg]", 0.06),
        ("theta_open_deg", "theta_closed_deg", "Pitch theta [deg]", 0.06),
        ("psi_open_deg", "psi_closed_deg", "Yaw psi [deg]", 0.5),
    ]

    for axis, (open_name, closed_name, ylabel, limit) in zip(
        axes,
        pairs,
    ):
        axis.plot(
            t_hr,
            data[open_name],
            label="Open loop",
        )

        axis.plot(
            t_hr,
            data[closed_name],
            "--",
            label="Closed loop",
        )

        axis.axhline(limit, linestyle=":")
        axis.axhline(-limit, linestyle=":")
        axis.set_ylabel(ylabel)
        axis.grid(True)
        axis.legend()

    axes[-1].set_xlabel("Time [hr]")

    fig.suptitle(
        f"HW4 {prefix}: Periodic Disturbance Response"
    )

    fig.tight_layout()

    # Control effort.
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(9, 8),
        sharex=True,
    )

    effort_data = [
        ("Tcx_Nm", "Tcx [N m]"),
        ("Tcy_Nm", "Tcy [N m]"),
        ("Tcz_Nm", "Tcz [N m]"),
    ]

    for axis, (name, ylabel) in zip(
        axes,
        effort_data,
    ):
        axis.plot(
            t_hr,
            data[name],
        )

        axis.set_ylabel(ylabel)
        axis.grid(True)

    axes[-1].set_xlabel("Time [hr]")

    fig.suptitle(
        f"HW4 {prefix}: Periodic-Disturbance Control Effort"
    )

    fig.tight_layout()


def plot_angle_commands(prefix):
    pitch = load_csv(
        f"hw4_{prefix}_pitch_command.csv"
    )

    roll = load_csv(
        f"hw4_{prefix}_roll_command.csv"
    )

    # Pitch command.
    fig, axis = plt.subplots(
        figsize=(9, 5),
    )

    axis.plot(
        pitch["time_s"],
        pitch["reference_deg"],
        "--",
        label="Reference",
    )

    axis.plot(
        pitch["time_s"],
        pitch["response_deg"],
        label="Response",
    )

    axis.axhline(0.06, linestyle=":")
    axis.set_xlabel("Time [s]")
    axis.set_ylabel("Pitch theta [deg]")
    axis.set_title(
        f"HW4 {prefix}: Pitch Angle Command"
    )
    axis.grid(True)
    axis.legend()

    fig.tight_layout()

    # Roll command + yaw cross coupling.
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
    )

    axes[0].plot(
        roll["time_s"],
        roll["reference_deg"],
        "--",
        label="Reference",
    )

    axes[0].plot(
        roll["time_s"],
        roll["response_deg"],
        label="Roll response",
    )

    axes[0].axhline(0.06, linestyle=":")
    axes[0].set_ylabel("Roll phi [deg]")
    axes[0].grid(True)
    axes[0].legend()

    # The C++ roll-command CSV does not store yaw cross response in
    # the same columns, so plot the effort channels here and leave
    # the yaw-cross curve to a future expanded CSV if desired.
    axes[1].plot(
        roll["time_s"],
        roll["effort1_Nm"],
        label="Tcx",
    )

    axes[1].plot(
        roll["time_s"],
        roll["effort2_Nm"],
        label="Tcz",
    )

    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Control torque [N m]")
    axes[1].grid(True)
    axes[1].legend()

    fig.suptitle(
        f"HW4 {prefix}: Roll Angle Command and Control Effort"
    )

    fig.tight_layout()


def plot_pitch_bode(prefix):
    data = load_csv(
        f"hw4_{prefix}_pitch_bode.csv"
    )

    frequency = data["omega_rad_s"]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
    )

    axes[0].semilogx(
        frequency,
        20.0 * np.log10(data["magnitude"]),
    )

    axes[0].axhline(0.0, linestyle=":")
    axes[0].set_ylabel("Magnitude [dB]")
    axes[0].grid(True, which="both")

    axes[1].semilogx(
        frequency,
        data["phase_deg"],
    )

    axes[1].axhline(-180.0, linestyle=":")
    axes[1].set_xlabel("Frequency [rad/s]")
    axes[1].set_ylabel("Phase [deg]")
    axes[1].grid(True, which="both")

    fig.suptitle(
        f"HW4 {prefix}: Pitch Open-Loop Frequency Response"
    )

    fig.tight_layout()


def plot_root_loci(prefix):
    pitch = load_csv(
        f"hw4_{prefix}_pitch_root_locus.csv"
    )

    ry = load_csv(
        f"hw4_{prefix}_rollyaw_root_locus.csv"
    )

    fig, axis = plt.subplots(
        figsize=(8, 7),
    )

    axis.plot(
        pitch["real1"],
        pitch["imag1"],
    )

    axis.plot(
        pitch["real2"],
        pitch["imag2"],
    )

    axis.axvline(0.0, linestyle=":")
    axis.axhline(0.0, linestyle=":")
    axis.set_xlabel("Real axis [s^-1]")
    axis.set_ylabel("Imaginary axis [s^-1]")
    axis.set_title(
        f"HW4 {prefix}: Pitch Root Locus"
    )
    axis.grid(True)

    fig.tight_layout()

    fig, axis = plt.subplots(
        figsize=(8, 7),
    )

    for suffix in [
        "1", "2", "3", "4"
    ]:
        axis.plot(
            ry[f"root{suffix}_real"],
            ry[f"root{suffix}_imag"],
        )

    axis.axvline(0.0, linestyle=":")
    axis.axhline(0.0, linestyle=":")
    axis.set_xlabel("Real axis [s^-1]")
    axis.set_ylabel("Imaginary axis [s^-1]")
    axis.set_title(
        f"HW4 {prefix}: Roll-Yaw Root Locus / Pole Sweep"
    )
    axis.grid(True)

    fig.tight_layout()


def main():
    variant = "design_a"

    if len(sys.argv) > 1:
        variant = sys.argv[1].lower()

    if variant not in {"design_a", "design_b"}:
        raise SystemExit(
            "Usage: python3 validation/hw4/visualize.py [design_a|design_b]"
        )

    plot_step_response(variant)
    plot_periodic(variant)
    plot_angle_commands(variant)
    plot_pitch_bode(variant)
    plot_root_loci(variant)

    plt.show()


if __name__ == "__main__":
    main()
