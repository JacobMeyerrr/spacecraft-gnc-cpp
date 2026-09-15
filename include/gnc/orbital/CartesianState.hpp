#pragma once

#include "gnc/math/Vector.hpp"

namespace gnc::orbital {

    struct CartesianState
    {
        gnc::math::Vector3 position;
        gnc::math::Vector3 velocity;
        double time;
    };

}