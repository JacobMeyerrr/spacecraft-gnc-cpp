#include "gnc/math/Polynomial.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gnc::math {

    std::vector<std::complex<double>> polynomialRoots(
        const std::vector<double>& coefficients)
    {
        if(coefficients.size() < 2 ||
           std::abs(coefficients.front()) < 1.0e-30)
        {
            throw std::invalid_argument(
                "Polynomial must have a nonzero leading coefficient"
            );
        }

        const int degree =
            static_cast<int>(coefficients.size()) - 1;

        Eigen::MatrixXcd companion =
            Eigen::MatrixXcd::Zero(degree, degree);

        for(int i = 1; i < degree; ++i)
        {
            companion(i, i - 1) = 1.0;
        }

        for(int i = 0; i < degree; ++i)
        {
            companion(i, degree - 1) =
                -coefficients[degree - i]
                / coefficients.front();
        }

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(
            companion
        );

        std::vector<std::complex<double>> roots;
        roots.reserve(degree);

        for(int i = 0; i < degree; ++i)
        {
            roots.push_back(
                solver.eigenvalues()(i)
            );
        }

        return roots;
    }

    double maxAbs(
        const std::vector<double>& values)
    {
        double result = 0.0;

        for(const double value : values)
        {
            result = std::max(
                result,
                std::abs(value)
            );
        }

        return result;
    }

}
