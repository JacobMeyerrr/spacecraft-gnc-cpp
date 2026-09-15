// Shared spacecraft GNC constants.
//
// The project uses SI units unless a validation case explicitly
// states otherwise.
#pragma once

namespace gnc::constants {

    // Mathematical constants
    inline constexpr double pi = 3.141592653589793238462643383279502884;
    inline constexpr double two_pi = 2.0 * pi;

    // Angle conversions
    inline constexpr double d2r = pi / 180.0;
    inline constexpr double r2d = 180.0 / pi;

    // Earth parameters
    inline constexpr double r_E = 6'378'137.0;      // m
    inline constexpr double mu_E = 3.986004418e14;  // m^3/s^2
    inline constexpr double n_E = 7.2921150e-5;     // rad/s
    inline constexpr double g0 = 9.80665;           // m/s^2
    inline constexpr double sidereal_day = 86'164.0905; // s

    // Numerical tolerance used by reusable numerical routines.
    inline constexpr double num_tolerance = 1.0e-12;

}
