#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ================================================================
# HW2 Problem 2 plotting
#
# C++ does the actual orbit/attitude calculation.
# Python is only used to visualize the resulting CSV.
# ================================================================


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_FILE = (
    PROJECT_ROOT
    / "data"
    / "validation"
    / "hw2"
    / "problem2_output.csv"
)


def main():

    if not DATA_FILE.exists():
        raise FileNotFoundError(
            "Run ./build/hw2_attitude first."
        )

    data = np.genfromtxt(
        DATA_FILE,
        delimiter=",",
        names=True
    )

    orbit_number = data["orbit_number"]


    # ============================================================
    # Figure 1: Euler angles
    # ============================================================

    roll = data["roll_deg"]
    pitch = data["pitch_deg"]
    yaw = data["yaw_deg"]

    roll_unwrapped = np.degrees(
        np.unwrap(
            np.radians(roll)
        )
    )

    yaw_unwrapped = np.degrees(
        np.unwrap(
            np.radians(yaw)
        )
    )

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(9, 8),
        sharex=True
    )

    axes[0].plot(
        orbit_number,
        yaw_unwrapped
    )

    axes[0].set_ylabel(
        "Yaw [deg]"
    )

    axes[0].grid(True)

    axes[1].plot(
        orbit_number,
        pitch
    )

    axes[1].set_ylabel(
        "Pitch [deg]"
    )

    axes[1].grid(True)

    axes[2].plot(
        orbit_number,
        roll_unwrapped
    )

    axes[2].set_ylabel(
        "Roll [deg]"
    )

    axes[2].set_xlabel(
        "Time [orbital periods]"
    )

    axes[2].grid(True)

    figure.suptitle(
        "Satellite 3-2-1 Euler Angles Relative to ECI"
    )

    figure.tight_layout()


    # ============================================================
    # Figure 2: Quaternion components
    # ============================================================

    figure, axis = plt.subplots(
        figsize=(9, 5)
    )

    axis.plot(
        orbit_number,
        data["q1"],
        label="q1"
    )

    axis.plot(
        orbit_number,
        data["q2"],
        label="q2"
    )

    axis.plot(
        orbit_number,
        data["q3"],
        label="q3"
    )

    axis.plot(
        orbit_number,
        data["q4"],
        label="q4"
    )

    axis.set_xlabel(
        "Time [orbital periods]"
    )

    axis.set_ylabel(
        "Quaternion component"
    )

    axis.set_title(
        "Body-to-ECI Quaternion Components"
    )

    axis.grid(True)
    axis.legend()

    figure.tight_layout()


    # ============================================================
    # Figure 3: Body axes as ECI components
    #
    # One figure per body axis, each with three components.
    # ============================================================

    body_axes = {
        "x_B": (
            data["body_x_eci_1"],
            data["body_x_eci_2"],
            data["body_x_eci_3"],
        ),
        "y_B": (
            data["body_y_eci_1"],
            data["body_y_eci_2"],
            data["body_y_eci_3"],
        ),
        "z_B": (
            data["body_z_eci_1"],
            data["body_z_eci_2"],
            data["body_z_eci_3"],
        ),
    }


    for axis_name, components in body_axes.items():

        figure, axes = plt.subplots(
            3,
            1,
            figsize=(9, 7),
            sharex=True
        )

        for component in range(3):

            axes[component].plot(
                orbit_number,
                components[component]
            )

            axes[component].set_ylabel(
                f"{axis_name} / ECI {component + 1}"
            )

            axes[component].grid(True)

        axes[-1].set_xlabel(
            "Time [orbital periods]"
        )

        figure.suptitle(
            f"{axis_name} Expressed in ECI"
        )

        figure.tight_layout()


    # ============================================================
    # Figure 4: 3D orbit + body axes
    #
    # Only one revolution is shown geometrically because the
    # ideal two-body orbit repeats every period.
    # ============================================================

    first_orbit = orbit_number <= 1.0

    x = data["x_eci_km"]
    y = data["y_eci_km"]
    z = data["z_eci_km"]


    figure = plt.figure(
        figsize=(10, 8)
    )

    axis = figure.add_subplot(
        111,
        projection="3d"
    )

    axis.plot(
        x[first_orbit],
        y[first_orbit],
        z[first_orbit]
    )


    indices = np.where(first_orbit)[0]

    snapshot_indices = np.linspace(
        indices[0],
        indices[-1],
        12,
        dtype=int
    )


    orbit_radius_scale = np.mean(
        np.sqrt(
            x[first_orbit] ** 2
            + y[first_orbit] ** 2
            + z[first_orbit] ** 2
        )
    )

    arrow_scale = 0.14 * orbit_radius_scale


    for index in snapshot_indices:

        origin = np.array([
            x[index],
            y[index],
            z[index]
        ])


        body_x = np.array([
            data["body_x_eci_1"][index],
            data["body_x_eci_2"][index],
            data["body_x_eci_3"][index]
        ])

        body_y = np.array([
            data["body_y_eci_1"][index],
            data["body_y_eci_2"][index],
            data["body_y_eci_3"][index]
        ])

        body_z = np.array([
            data["body_z_eci_1"][index],
            data["body_z_eci_2"][index],
            data["body_z_eci_3"][index]
        ])


        axis.quiver(
            origin[0],
            origin[1],
            origin[2],
            *(arrow_scale * body_x)
        )

        axis.quiver(
            origin[0],
            origin[1],
            origin[2],
            *(arrow_scale * body_y)
        )

        axis.quiver(
            origin[0],
            origin[1],
            origin[2],
            *(arrow_scale * body_z)
        )


    axis.set_xlabel(
        "ECI X [km]"
    )

    axis.set_ylabel(
        "ECI Y [km]"
    )

    axis.set_zlabel(
        "ECI Z [km]"
    )

    axis.set_title(
        "Earth-Pointing Spacecraft Orbit and Body Attitude"
    )

    axis.grid(True)

    figure.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()