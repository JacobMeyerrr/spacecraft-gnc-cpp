# Build and Run

From the repository root:

```bash
rm -rf build
cmake -S . -B build
cmake --build build

./build/validate_hw1_orbit
./build/validate_hw1_attitude
./build/validate_hw2_quaternions
./build/validate_hw2_attitude
./build/validate_hw2_pointing
./build/validate_hw3
./build/validate_hw4
```

Python visualizations:

```bash
python3 validation/hw1/visualize_problem1.py
python3 validation/hw2/visualize_problem2.py
python3 validation/hw2/visualize_problem3.py
python3 validation/hw3/visualize.py
python3 validation/hw4/visualize.py
```

The reusable library is the `spacecraft_gnc` C++ library. HW1-HW4 are validation clients, not the architecture of the library.
