#include "gnc/dynamics/RigidBodyDynamics.hpp"

namespace gnc::dynamics {

    gnc::math::Vector3 angularAcceleration(
        const gnc::math::Vector3& omega_body,
        const gnc::math::Vector3& torque_body,
        const gnc::math::Matrix3& inertia)
    {
        return inertia.inverse() * (
            torque_body
            - omega_body.cross(inertia * omega_body)
        );
    }

}
