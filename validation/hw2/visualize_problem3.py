#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ================================================================
# HW2 Problem 3 visualization
# ================================================================


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_FILE = (
    PROJECT_ROOT
    / "data"
    / "validation"
    / "hw2"
    / "problem3_output.csv"
)


def main():

    if not DATA_FILE.exists():
        raise FileNotFoundError(
            "Run ./build/hw2_pointing first."
        )

    data = np.genfromtxt(
        DATA_FILE,
        delimiter=",",
        names=True
    )

    orbit_number = data["orbit_number"]


    # ============================================================
    # Problem 3(a): commanded Euler angles
    # ============================================================

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(9, 8),
        sharex=True
    )

    axes[0].plot(
        orbit_number,
        data["roll_deg"]
    )

    axes[0].set_ylabel(
        "Roll phi [deg]"
    )

    axes[0].grid(True)


    axes[1].plot(
        orbit_number,
        data["pitch_deg"]
    )

    axes[1].set_ylabel(
        "Pitch theta [deg]"
    )

    axes[1].grid(True)


    axes[2].plot(
        orbit_number,
        data["yaw_deg"]
    )

    axes[2].set_ylabel(
        "Yaw psi [deg]"
    )

    axes[2].set_xlabel(
        "Time [orbital periods]"
    )

    axes[2].grid(True)


    figure.suptitle(
        "XPOP-YPSL Commanded 3-2-1 Euler Angles"
    )

    figure.tight_layout()


    # ============================================================
    # Problem 3(b): solar array gimbal
    # ============================================================

    figure, axis = plt.subplots(
        figsize=(9, 5)
    )

    axis.plot(
        orbit_number,
        data["delta_deg"]
    )

    axis.set_xlabel(
        "Time [orbital periods]"
    )

    axis.set_ylabel(
        "Solar array gimbal delta [deg]"
    )

    axis.set_title(
        "Solar Array Gimbal Angle"
    )

    axis.grid(True)

    figure.tight_layout()


    # ============================================================
    # Problem 3(c): extra-credit magnetic field direction
    # ============================================================

    figure, axis = plt.subplots(
        figsize=(9, 5)
    )

    axis.plot(
        orbit_number,
        data["Bx_body"],
        label="B_xB"
    )

    axis.plot(
        orbit_number,
        data["By_body"],
        label="B_yB"
    )

    axis.plot(
        orbit_number,
        data["Bz_body"],
        label="B_zB"
    )

    axis.set_xlabel(
        "Time [orbital periods]"
    )

    axis.set_ylabel(
        "Unit magnetic-field component"
    )

    axis.set_ylim(
        -1.0,
        1.0
    )

    axis.set_title(
        "Unit Magnetic Field Direction Observed in Body Frame"
    )

    axis.grid(True)
    axis.legend()

    figure.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()