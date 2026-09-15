#include "gnc/control/StateSpace.hpp"

#include <Eigen/Eigenvalues>

namespace gnc::control {

    std::vector<std::complex<double>> StateSpaceModel::poles() const
    {
        Eigen::EigenSolver<gnc::math::Matrix4> solver(A);
        std::vector<std::complex<double>> result;
        result.reserve(4);

        for(int i = 0; i < 4; ++i)
        {
            result.push_back(solver.eigenvalues()(i));
        }

        return result;
    }

}
