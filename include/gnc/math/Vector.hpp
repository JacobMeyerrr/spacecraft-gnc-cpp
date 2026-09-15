// Vector aliases used throughout the GNC library.
#pragma once

#include <Eigen/Dense>

namespace gnc::math {

    using Vector2 = Eigen::Vector2d;
    using Vector3 = Eigen::Vector3d;
    using Vector4 = Eigen::Vector4d;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

}
