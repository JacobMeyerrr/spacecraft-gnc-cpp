#pragma once

#include "gnc/math/Matrix.hpp"

namespace gnc::dynamics {

    struct SpacecraftProperties
    {
        gnc::math::Matrix3 inertia = gnc::math::Matrix3::Identity();

        gnc::math::Matrix3 inverseInertia() const;
    };

}
