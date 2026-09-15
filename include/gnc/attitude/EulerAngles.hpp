#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::attitude {

    struct EulerAngles
    {
        double phi{0.0};
        double theta{0.0};
        double psi{0.0};
    };

    // Course 1-2-3 body-to-reference DCM.
    gnc::math::Matrix3 dcm123(
        const EulerAngles& angles);

    // Course 3-2-1 Euler-angle extraction from the active DCM
    // convention used in the reusable attitude library.
    EulerAngles dcmToEuler321(
        const gnc::math::Matrix3& R);

    // Construct the equivalent 3-2-1 DCM.
    gnc::math::Matrix3 dcm321(
        const EulerAngles& angles);

}
