# Design Notes

## Purpose

The repository is a reusable spacecraft GNC library. Validation cases are kept outside the core library so that adding an assignment, textbook case, mission scenario, or regression test does not force the reusable source tree to follow that case.

## Package responsibilities

### attitude

Rotation matrices, Euler angle parameterizations, differential kinematics, and quaternions.

### orbital

Classical orbital elements, two-body orbit calculations, inertial propagation, and reusable orbital reference-frame construction.

### dynamics

Rigid-body rotational dynamics, inertia handling, and gravity-gradient torque.

### models

Composable spacecraft-level models that assemble lower-level attitude, orbital, and dynamics functionality.

### control

PID/PD primitives, state-space models, characteristic polynomials, and momentum-biased attitude-controller design.

### simulation

Numerical integration, reusable fixed-step simulation, and CSV logging.

## Validation layers

HW1 checks the basic orbit and rotation primitives.

HW2 checks quaternion operations, Earth-pointing reference frames, and XPOP/-YPSL attitude calculations.

HW3 checks nonlinear spacecraft rotational dynamics against its linearized perturbation model and transfer-function behavior.

HW4 checks momentum-biased controller design, coupled roll/yaw dynamics, pitch control, disturbance rejection, angle commands, and frequency-domain margins.

The validation code should call the library rather than duplicate the implementation of the physical models it is intended to validate.
