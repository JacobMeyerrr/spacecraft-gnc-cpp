#pragma once

#include <complex>
#include <vector>

namespace gnc::math {

    // Return the roots of a polynomial whose coefficients are ordered:
    //
    //     a_n, a_(n-1), ..., a_1, a_0
    //
    // using a companion-matrix eigenvalue calculation.
    std::vector<std::complex<double>> polynomialRoots(
        const std::vector<double>& coefficients);

    double maxAbs(
        const std::vector<double>& values);

}
