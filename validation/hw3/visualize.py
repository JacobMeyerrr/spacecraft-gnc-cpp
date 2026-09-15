#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HW3 = ROOT / "data" / "validation" / "hw3"


def load_case(case_id):
    path = HW3 / f"hw3_case{case_id}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Run ./build/validate_hw3 first: {path}")
    return np.genfromtxt(path, delimiter=",", names=True)


def plot_case(data, case_id):
    t = data["time_s"]
    names = ["phi", "theta", "psi"]
    labels = ["Roll phi [deg]", "Pitch theta [deg]", "Yaw psi [deg]"]

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for i, (name, label) in enumerate(zip(names, labels)):
        axes[i].plot(t, data[f"{name}_linear_deg"], label="Linear")
        axes[i].plot(t, data[f"{name}_nonlinear_deg"], "--", label="Nonlinear")
        axes[i].set_ylabel(label)
        axes[i].grid(True)
        if i == 0:
            axes[i].legend()
    axes[-1].set_xlabel("Time [s]")
    fig.suptitle(f"HW3 Case {case_id}: Linear vs Nonlinear Attitude")
    fig.tight_layout()

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for i, (name, label) in enumerate(zip(names, labels)):
        axes[i].plot(t, data[f"{name}_error_deg"])
        axes[i].set_ylabel(f"{name} error [deg]")
        axes[i].grid(True)
    axes[-1].set_xlabel("Time [s]")
    fig.suptitle(f"HW3 Case {case_id}: Nonlinear - Linear")
    fig.tight_layout()

    if case_id == 3:
        fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
        torque_labels = ["Tx [N m]", "Ty [N m]", "Tz [N m]"]
        torque_names = ["Tx_Nm", "Ty_Nm", "Tz_Nm"]
        for i, (name, label) in enumerate(zip(torque_names, torque_labels)):
            axes[i].plot(t, data[name])
            axes[i].set_ylabel(label)
            axes[i].grid(True)
        axes[-1].set_xlabel("Time [s]")
        fig.suptitle("HW3 Case 3: Input Disturbance Torque")
        fig.tight_layout()


def plot_frequency_response():
    path = HW3 / "hw3_frequency_response.csv"
    if not path.exists():
        return
    data = np.genfromtxt(path, delimiter=",", names=True)

    names = [
        "phi_Tx_mag", "phi_Ty_mag", "phi_Tz_mag",
        "theta_Tx_mag", "theta_Ty_mag", "theta_Tz_mag",
        "psi_Tx_mag", "psi_Ty_mag", "psi_Tz_mag",
    ]
    labels = [
        "phi/Tx", "phi/Ty", "phi/Tz",
        "theta/Tx", "theta/Ty", "theta/Tz",
        "psi/Tx", "psi/Ty", "psi/Tz",
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    for name, label in zip(names, labels):
        ax.loglog(data["omega_rad_s"], data[name], label=label)
    ax.set_xlabel("Frequency omega [rad/s]")
    ax.set_ylabel("Transfer magnitude")
    ax.set_title("HW3 Torque-to-Attitude Transfer Matrix Magnitudes")
    ax.grid(True, which="both")
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()


def main():
    for case_id in (1, 2, 3):
        plot_case(load_case(case_id), case_id)
    plot_frequency_response()
    plt.show()


if __name__ == "__main__":
    main()
