HW4 periodic-validation fix
===========================

The HW4 periodic roll-yaw validation path was rewritten to integrate the
four scalar states [phi, phi_dot, psi, psi_dot] explicitly with RK4 using
the Lecture 18 equations. The previous path relied on an Eigen fixed-size
state inside the periodic validation loop and produced invalid near-zero
outputs on the local run.

Expected Design B regression values after rebuilding/running:

Closed-loop periodic attitude maxima:
  roll  ~= 0.05087 deg
  pitch ~= 0.03000 deg
  yaw   ~= 0.39436 deg

Periodic control-effort maxima:
  |Tcx| ~= 1.794e-2 N m
  |Tcy| ~= 1.362e-4 N m
  |Tcz| ~= 4.416e-4 N m

The design parameters remain:
  HB = 1200 N m s
  K  = 20 N m/rad
  tau = 22.36067977 s
  k  = 0.0246174146

The assignment specifies the simultaneous periodic disturbances and the
three axis accuracy requirements; the reference design uses the Lecture 18
roll-yaw controller. The corresponding expected results are documented in
HW4 solutions included with the project.
