#pragma once

#include "gnc/orbital/CartesianState.hpp"
#include "gnc/orbital/Orbit.hpp"

namespace gnc::orbital {

    CartesianState propagateTwoBody(
        const Orbit& orbit,
        double time);

}
