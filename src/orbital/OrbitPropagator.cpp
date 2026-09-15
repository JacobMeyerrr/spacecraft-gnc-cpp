#include "gnc/math/Matrix.hpp"
#include "gnc/orbital/OrbitPropagator.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Vector.hpp"

#include <cmath>

namespace gnc::orbital {

    CartesianState propagateTwoBody(
        const Orbit& orbit,
        double time)
    {
        const auto& e = orbit.elements();

        const double a = e.semi_major_axis;
        const double ecc = e.eccentricity;
        const double i = e.inclination;
        const double raan = e.raan;
        const double argp = e.argument_of_periapsis;
        const double nu = orbit.trueAnomaly(time);

        const double p =
            a * (1.0 - ecc * ecc);

        const double radius =
            p / (1.0 + ecc * std::cos(nu));

        gnc::math::Vector3 r_pqw;
        r_pqw <<
            radius * std::cos(nu),
            radius * std::sin(nu),
            0.0;

        gnc::math::Vector3 v_pqw;
        const double scale =
            std::sqrt(gnc::constants::mu_E / p);

        v_pqw <<
            -scale * std::sin(nu),
            scale * (ecc + std::cos(nu)),
            0.0;

        const auto R1 = [](double angle)
        {
            const double c = std::cos(angle);
            const double s = std::sin(angle);

            return (gnc::math::Matrix3() <<
                1.0, 0.0, 0.0,
                0.0, c, -s,
                0.0, s, c
            ).finished();
        };

        const auto R3 = [](double angle)
        {
            const double c = std::cos(angle);
            const double s = std::sin(angle);

            return (gnc::math::Matrix3() <<
                c, -s, 0.0,
                s, c, 0.0,
                0.0, 0.0, 1.0
            ).finished();
        };

        const gnc::math::Matrix3 Q =
            R3(raan) *
            R1(i) *
            R3(argp);

        CartesianState state;
        state.position = Q * r_pqw;
        state.velocity = Q * v_pqw;
        state.time = time;

        return state;
    }

}
