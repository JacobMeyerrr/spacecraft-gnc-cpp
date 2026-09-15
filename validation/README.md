# Validation

These programs are regression cases for the reusable spacecraft GNC library.

They intentionally live outside `include/gnc/` and `src/` so the core architecture remains independent of any single assignment.

- `hw1/` — orbit and Euler rotation checks
- `hw2/` — quaternion, Earth-pointing attitude, and XPOP/-YPSL checks
- `hw3/` — linear/nonlinear rigid-body dynamics checks
- `hw4/` — momentum-biased controller design and verification
