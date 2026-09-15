# Spacecraft GNC C++

A modular C++ spacecraft Guidance, Navigation, and Control codebase with a current emphasis on attitude determination and control.

The design separates reusable spacecraft/GNC functionality from assignment-specific validation. Homework 1 through Homework 4 live under `validation/` and act as regression cases for the library rather than defining its architecture.

## Architecture

```text
spacecraft-gnc-cpp/
├── include/gnc/
│   ├── attitude/       # DCMs, Euler angles, quaternions, kinematics
│   ├── control/        # PID/PD, state-space, controller design
│   ├── dynamics/       # rigid-body and environmental dynamics
│   ├── math/           # Eigen aliases, constants, numerics, polynomials
│   ├── models/         # reusable spacecraft-level models
│   ├── orbital/        # orbital mechanics and reference frames
│   └── simulation/     # integration, simulation, logging
│
├── src/                # Implementations of the reusable library
│
├── validation/
│   ├── hw1/            # Orbit + Euler rotation regression cases
│   ├── hw2/            # Quaternion + Earth-pointing + XPOP/YPSL cases
│   ├── hw3/            # Linear/nonlinear spacecraft dynamics cases
│   └── hw4/            # Momentum-biased control design cases
│
├── data/validation/    # Generated CSV outputs from validation runs
└── docs/               # Design notes and architectural documentation
```

## Units and conventions

The reusable library uses SI units by default. The HW2 orbital validation keeps its historical kilometer-based constants inside that validation case because the original numerical reference uses kilometers.

The quaternion library uses the course-compatible vector-first, scalar-last representation:

```text
q = [q1 q2 q3 q4]^T
```

with `q1:q3` as the vector part and `q4` as the scalar part.

## Build

Eigen3 is required.

```bash
cmake -S . -B build
cmake --build build
```

## Run validations

```bash
./build/validate_hw1_orbit
./build/validate_hw1_attitude

./build/validate_hw2_quaternions
./build/validate_hw2_attitude
./build/validate_hw2_pointing

./build/validate_hw3
./build/validate_hw4
```

HW4 accepts `design_a`, `design_b`, or `all`:

```bash
./build/validate_hw4 design_b
./build/validate_hw4 design_a
./build/validate_hw4 all
```

## Visualization

```bash
python3 validation/hw1/visualize_problem1.py
python3 validation/hw2/visualize_problem2.py
python3 validation/hw2/visualize_problem3.py
python3 validation/hw3/visualize.py
python3 validation/hw4/visualize.py design_b
```

## Validation philosophy

The homework cases are useful regression tests because they have known expected outputs, but they are not the final proof of correctness of the reusable library. Independent textbook/reference cases should be added under a separate validation area before treating the attitude and control conventions as fully certified.
