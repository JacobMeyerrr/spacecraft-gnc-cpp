// Matrix aliases used throughout the GNC library.
#pragma once

#include <Eigen/Dense>

namespace gnc::math {

    using Matrix2 = Eigen::Matrix2d;
    using Matrix3 = Eigen::Matrix3d;
    using Matrix4 = Eigen::Matrix4d;
    using Matrix6 = Eigen::Matrix<double, 6, 6>;

}
