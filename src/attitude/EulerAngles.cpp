#include "gnc/attitude/EulerAngles.hpp"
#include "gnc/math/Numerics.hpp"

#include <cmath>

namespace gnc::attitude {

    gnc::math::Matrix3 dcm123(
        const EulerAngles& a)
    {
        const double cphi = std::cos(a.phi);
        const double sphi = std::sin(a.phi);
        const double ctheta = std::cos(a.theta);
        const double stheta = std::sin(a.theta);
        const double cpsi = std::cos(a.psi);
        const double spsi = std::sin(a.psi);

        return (gnc::math::Matrix3() <<
            cpsi * ctheta,
            spsi * ctheta,
            -stheta,

            -cphi * spsi + sphi * stheta * cpsi,
            sphi * stheta * spsi + cphi * cpsi,
            sphi * ctheta,

            cphi * stheta * cpsi + sphi * spsi,
            cphi * stheta * spsi - sphi * cpsi,
            cphi * ctheta
        ).finished();
    }

    EulerAngles dcmToEuler321(
        const gnc::math::Matrix3& R)
    {
        EulerAngles angles;

        angles.theta = std::asin(
            gnc::math::clampUnit(-R(0,2))
        );

        angles.psi = std::atan2(
            R(0,1),
            R(0,0)
        );

        angles.phi = -std::atan2(
            R(1,2),
            R(2,2)
        );

        return angles;
    }

    gnc::math::Matrix3 dcm321(
        const EulerAngles& a)
    {
        const double cphi = std::cos(a.phi);
        const double sphi = std::sin(a.phi);
        const double ctheta = std::cos(a.theta);
        const double stheta = std::sin(a.theta);
        const double cpsi = std::cos(a.psi);
        const double spsi = std::sin(a.psi);

        return (gnc::math::Matrix3() <<
            cpsi * ctheta,
            spsi * ctheta,
            -stheta,

            -cphi * spsi + sphi * stheta * cpsi,
            cphi * cpsi + sphi * stheta * spsi,
            sphi * ctheta,

            sphi * spsi + cphi * stheta * cpsi,
            -sphi * cpsi + cphi * stheta * spsi,
            cphi * ctheta
        ).finished();
    }

}
