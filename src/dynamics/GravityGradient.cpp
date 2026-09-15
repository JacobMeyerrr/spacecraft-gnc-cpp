#include "gnc/dynamics/GravityGradient.hpp"

namespace gnc::dynamics {

    gnc::math::Vector3 gravityGradientTorque(
        double orbital_rate,
        const gnc::math::Vector3& nadir_body,
        const gnc::math::Matrix3& inertia)
    {
        return 3.0 * orbital_rate * orbital_rate
            * nadir_body.cross(inertia * nadir_body);
    }

}
