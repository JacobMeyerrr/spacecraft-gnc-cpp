#!/usr/bin/env python3

import math
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = PROJECT_ROOT / "build" / "validate_hw1_orbit"

# Professor's solution-case parameters.
# These are used for the extra diagnostic plots that are not printed by the C++ executable.
HP = 500_000.0
E = 0.9
INC = math.radians(30.0)
RAAN0 = math.radians(45.0)
ARG_PERIGEE0 = math.radians(-30.0)

# Keep these consistent with the professor's supplied solution for verification.
MU = 4.035040e14
RE = 6.3753845e6
J = 1.624e-3

RP = HP + RE
RA = RP * (1.0 + E) / (1.0 - E)
A = (RA + RP) / 2.0
N = math.sqrt(MU / A**3)
PERIOD = 2.0 * math.pi / N

A_RE_SQ = (A / RE)**2
ONE_MINUS_E2 = 1.0 - E**2
RAAN_RATE = -(J * N * math.cos(INC)) / (A_RE_SQ * ONE_MINUS_E2**2)
ARG_PERIGEE_RATE = (
    J * N * (4.0 - 5.0 * math.sin(INC)**2)
    / (2.0 * A_RE_SQ * ONE_MINUS_E2**2)
)


def R1(angle):
    s = math.sin(angle)
    c = math.cos(angle)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0, c, s],
        [0.0, -s, c],
    ])


def R2(angle):
    s = math.sin(angle)
    c = math.cos(angle)
    return np.array([
        [c, 0.0, -s],
        [0.0, 1.0, 0.0],
        [s, 0.0, c],
    ])


def R3(angle):
    s = math.sin(angle)
    c = math.cos(angle)
    return np.array([
        [c, s, 0.0],
        [-s, c, 0.0],
        [0.0, 0.0, 1.0],
    ])


def run_cpp():
    if not EXECUTABLE.exists():
        raise FileNotFoundError(
            f"{EXECUTABLE} does not exist. Build hw1_orbit first."
        )

    result = subprocess.run(
        [str(EXECUTABLE)],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def parse_cpp_output(output):
    rows = []
    reading = False

    for line in output.splitlines():
        line = line.strip()

        if line == "Time [s], True Anomaly [rad], Velocity [m/s]":
            reading = True
            continue

        if not reading or not line:
            continue

        parts = [x.strip() for x in line.split(",")]

        if len(parts) != 3:
            continue

        try:
            rows.append(tuple(float(x) for x in parts))
        except ValueError:
            pass

    if not rows:
        raise ValueError("Could not parse orbital data from hw1_orbit.")

    return np.array(rows)


def solve_kepler(M, eccentricity, guess=0.0):
    E_anomaly = guess

    for _ in range(100):
        f = E_anomaly - eccentricity * math.sin(E_anomaly) - M
        df = 1.0 - eccentricity * math.cos(E_anomaly)
        delta = f / df
        E_anomaly -= delta

        if abs(delta) < 1e-12:
            break

    return E_anomaly


def quaternion_from_matrix(R):
    # Vector-first, scalar-last convention: [q1, q2, q3, q4].
    q4_sq = max(0.0, 0.25 * (1.0 + np.trace(R)))
    q4 = math.sqrt(q4_sq)

    if abs(q4) > 1e-8:
        q1 = (R[1, 2] - R[2, 1]) / (4.0 * q4)
        q2 = (R[2, 0] - R[0, 2]) / (4.0 * q4)
        q3 = (R[0, 1] - R[1, 0]) / (4.0 * q4)
    else:
        # Eigen-decomposition fallback near a 180-degree rotation.
        eigvals, eigvecs = np.linalg.eig(R)
        axis = np.real(eigvecs[:, np.argmin(np.abs(eigvals - 1.0))])
        axis /= np.linalg.norm(axis)
        q1, q2, q3 = axis

    q = np.array([q1, q2, q3, q4])
    norm = np.linalg.norm(q)

    if norm != 0.0:
        q /= norm

    return q


def axis_angle_from_matrix(R):
    value = np.clip(0.5 * (np.trace(R) - 1.0), -1.0, 1.0)
    angle = math.acos(value)
    s = math.sin(angle)

    if abs(s) < 1e-10:
        axis = np.zeros(3)
    else:
        axis = np.array([
            R[1, 2] - R[2, 1],
            R[2, 0] - R[0, 2],
            R[0, 1] - R[1, 0],
        ]) / (2.0 * s)

    return angle, axis


def solution_case_history():
    # Match the professor's MATLAB solution: T/50, including both endpoints.
    times = np.linspace(0.0, PERIOD, 51)

    eccentric_anomaly = []
    true_anomaly = []
    radius = []
    altitude = []
    velocity = []
    angular_velocity = []

    euler313 = []
    quaternions = []
    rotation_angles = []
    rotation_axes = []

    previous_E = 0.0

    ROR = np.array([
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [-1.0, 0.0, 0.0],
    ])

    # Professor's supplied solution uses pitch alignment error = 15 deg here.
    phi = 0.0
    theta_alignment = math.radians(15.0)
    psi = 0.0
    RRB = R1(phi) @ R2(theta_alignment) @ R3(psi)

    for t in times:
        M = N * t
        E_now = solve_kepler(M, E, previous_E)
        previous_E = E_now

        r = A * (1.0 - E * math.cos(E_now))
        h = r - RE

        nu = 2.0 * math.atan2(
            math.sqrt(1.0 + E) * math.sin(E_now / 2.0),
            math.sqrt(1.0 - E) * math.cos(E_now / 2.0),
        )

        v = math.sqrt(MU * (2.0 / r - 1.0 / A))
        omega = math.sqrt(MU * A * (1.0 - E**2)) / r**2

        raan = RAAN0 + RAAN_RATE * t
        arg_perigee = ARG_PERIGEE0 + ARG_PERIGEE_RATE * t
        alpha = arg_perigee + nu

        RIN = R1(INC) @ R3(raan)
        RNO = R3(alpha)
        RIB = RIN @ RNO @ ROR @ RRB

        # Match professor's Euler-313 extraction formulas.
        p1 = -math.atan2(RIB[2, 0], RIB[2, 1])
        p2 = math.acos(np.clip(-RIB[2, 2], -1.0, 1.0))
        p3 = math.atan2(RIB[0, 2], RIB[1, 2])

        q = quaternion_from_matrix(RIB)
        rotation_angle, rotation_axis = axis_angle_from_matrix(RIB)

        eccentric_anomaly.append(E_now)
        true_anomaly.append(nu)
        radius.append(r)
        altitude.append(h)
        velocity.append(v)
        angular_velocity.append(omega)
        euler313.append([p1, p2, p3])
        quaternions.append(q)
        rotation_angles.append(rotation_angle)
        rotation_axes.append(rotation_axis)

    return {
        "time": times,
        "time_orbits": times / PERIOD,
        "E": np.array(eccentric_anomaly),
        "nu": np.array(true_anomaly),
        "r": np.array(radius),
        "h": np.array(altitude),
        "v": np.array(velocity),
        "omega": np.array(angular_velocity),
        "euler313": np.array(euler313),
        "q": np.array(quaternions),
        "axis_angle": np.array(rotation_angles),
        "axis": np.array(rotation_axes),
    }


def make_plots(data):
    t = data["time_orbits"]

    # 1. Spacecraft orbit
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    ax.plot(data["nu"], data["r"] / 1000.0, "o")
    ax.set_title("Spacecraft Orbit")
    fig.tight_layout()

    # 2. Eccentric anomaly
    plt.figure()
    plt.plot(t, np.degrees(data["E"]), linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("E [deg]")
    plt.title("Eccentric Anomaly E")
    plt.grid(True)
    plt.tight_layout()

    # 3. True anomaly
    plt.figure()
    plt.plot(t, np.degrees(data["nu"]), linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("True Anomaly [deg]")
    plt.title("True Anomaly")
    plt.grid(True)
    plt.tight_layout()

    # 4. Altitude
    plt.figure()
    plt.plot(t, data["h"] / 1000.0, linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("h [km]")
    plt.title("Spacecraft Altitude")
    plt.grid(True)
    plt.tight_layout()

    # 5. Linear velocity
    plt.figure()
    plt.plot(t, data["v"] / 1000.0, linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("v [km/s]")
    plt.title("Spacecraft Linear Velocity")
    plt.grid(True)
    plt.tight_layout()

    # 6. Angular velocity
    plt.figure()
    plt.plot(t, data["omega"], linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("omega [rad/s]")
    plt.title("Spacecraft Angular Velocity")
    plt.grid(True)
    plt.tight_layout()

    # 7. Euler 313
    plt.figure()
    euler_deg = np.degrees(data["euler313"])
    plt.plot(t, euler_deg[:, 0], label="Roll", linewidth=2)
    plt.plot(t, euler_deg[:, 1], label="Pitch", linewidth=2)
    plt.plot(t, euler_deg[:, 2], label="Yaw", linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("Euler Angles [deg]")
    plt.title("Spacecraft Orientation: Euler313 Parameterization")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # 8. Quaternion
    plt.figure()
    plt.plot(t, data["q"][:, 0], label="q1", linewidth=2)
    plt.plot(t, data["q"][:, 1], label="q2", linewidth=2)
    plt.plot(t, data["q"][:, 2], label="q3", linewidth=2)
    plt.plot(t, data["q"][:, 3], label="q4", linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("Quaternion Elements")
    plt.title("Spacecraft Orientation: Quaternion Parameterization")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # 9a. Axis-angle rotation angle
    plt.figure()
    plt.plot(t, np.degrees(data["axis_angle"]), linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("Rotation Angle [deg]")
    plt.title("Spacecraft Orientation: Axis/Angle Rotation Angle")
    plt.grid(True)
    plt.tight_layout()

    # 9b. Axis-angle rotation vector
    plt.figure()
    plt.plot(t, data["axis"][:, 0], label="n1", linewidth=2)
    plt.plot(t, data["axis"][:, 1], label="n2", linewidth=2)
    plt.plot(t, data["axis"][:, 2], label="n3", linewidth=2)
    plt.xlabel("t [orbit]")
    plt.ylabel("Rotation Vector")
    plt.title("Spacecraft Orientation: Axis/Angle Rotation Axis")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.show()


def main():
    try:
        cpp_output = run_cpp()
        cpp_rows = parse_cpp_output(cpp_output)

        print(f"Parsed {len(cpp_rows)} samples from C++ hw1_orbit.")
        print("Generating professor-solution-case diagnostic plots...")

        data = solution_case_history()
        make_plots(data)

    except (
        FileNotFoundError,
        ValueError,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
