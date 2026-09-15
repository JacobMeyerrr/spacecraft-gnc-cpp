#include "gnc/dynamics/SpacecraftProperties.hpp"

namespace gnc::dynamics {

    gnc::math::Matrix3 SpacecraftProperties::inverseInertia() const
    {
        return inertia.inverse();
    }

}
