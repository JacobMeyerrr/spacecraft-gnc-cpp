#pragma once

#include "gnc/math/Vector.hpp"
#include "gnc/math/Matrix.hpp"

namespace gnc::attitude {

    // Exact differential kinematics for the course 1-2-3 Euler
    // parameterization used by HW3.
    gnc::math::Matrix3 differentialKinematics123(
        const gnc::math::Vector3& euler_angles);

    gnc::math::Vector3 eulerRatesFromBodyRate123(
        const gnc::math::Vector3& euler_angles,
        const gnc::math::Vector3& body_rate_reference);

    gnc::math::Vector3 bodyRateFromEulerRates123(
        const gnc::math::Vector3& euler_angles,
        const gnc::math::Vector3& euler_rates);

}
