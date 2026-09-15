#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::dynamics {

    gnc::math::Vector3 gravityGradientTorque(
        double orbital_rate,
        const gnc::math::Vector3& nadir_body,
        const gnc::math::Matrix3& inertia);

}
