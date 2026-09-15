#include "gnc/attitude/AttitudeKinematics.hpp"

#include <cmath>
#include <stdexcept>

namespace gnc::attitude {

    gnc::math::Matrix3 differentialKinematics123(
        const gnc::math::Vector3& euler_angles)
    {
        const double theta = euler_angles(1);
        const double psi = euler_angles(2);

        const double ctheta = std::cos(theta);
        const double stheta = std::sin(theta);
        const double cpsi = std::cos(psi);
        const double spsi = std::sin(psi);

        if(std::abs(ctheta) < 1.0e-10)
        {
            throw std::runtime_error(
                "123 Euler representation is near gimbal lock"
            );
        }

        return (1.0 / ctheta) * (gnc::math::Matrix3() <<
            ctheta, spsi * stheta, cpsi * stheta,
            0.0, cpsi * ctheta, -spsi * ctheta,
            0.0, spsi, cpsi
        ).finished();
    }

    gnc::math::Vector3 eulerRatesFromBodyRate123(
        const gnc::math::Vector3& euler_angles,
        const gnc::math::Vector3& body_rate_reference)
    {
        return differentialKinematics123(euler_angles)
            * body_rate_reference;
    }

    gnc::math::Vector3 bodyRateFromEulerRates123(
        const gnc::math::Vector3& euler_angles,
        const gnc::math::Vector3& euler_rates)
    {
        return differentialKinematics123(euler_angles)
            .inverse()
            * euler_rates;
    }

}
