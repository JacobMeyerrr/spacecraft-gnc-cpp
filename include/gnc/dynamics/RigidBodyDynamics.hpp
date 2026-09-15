#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::dynamics {

    gnc::math::Vector3 angularAcceleration(
        const gnc::math::Vector3& omega_body,
        const gnc::math::Vector3& torque_body,
        const gnc::math::Matrix3& inertia);

}
