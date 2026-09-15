#pragma once

#include "gnc/math/Matrix.hpp"

#include <complex>
#include <vector>

namespace gnc::control {

    struct StateSpaceModel
    {
        gnc::math::Matrix4 A = gnc::math::Matrix4::Zero();
        Eigen::Matrix<double, 4, 2> B = Eigen::Matrix<double, 4, 2>::Zero();
        Eigen::Matrix<double, 2, 4> C = Eigen::Matrix<double, 2, 4>::Zero();
        Eigen::Matrix2d D = Eigen::Matrix2d::Zero();

        std::vector<std::complex<double>> poles() const;
    };

    struct PitchTransferFunctions
    {
        std::vector<double> open_loop_num;
        std::vector<double> open_loop_den;
        std::vector<double> disturbance_closed_num;
        std::vector<double> disturbance_closed_den;
        std::vector<double> command_closed_num;
        std::vector<double> command_closed_den;
    };

}
