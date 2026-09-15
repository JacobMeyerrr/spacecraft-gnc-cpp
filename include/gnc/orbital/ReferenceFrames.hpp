#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::orbital {

    struct ReferenceFrame
    {
        gnc::math::Vector3 x_eci;
        gnc::math::Vector3 y_eci;
        gnc::math::Vector3 z_eci;
        gnc::math::Matrix3 eci_to_reference;
    };

    // Earth-pointing/LVLH-style frame used by HW2 and the ADCS
    // examples:
    //
    //   +x = direction of motion
    //   +y = negative orbit normal
    //   +z = nadir
    ReferenceFrame earthPointingFrame(
        const gnc::math::Vector3& position_eci,
        const gnc::math::Vector3& velocity_eci);

}
