#include "gnc/orbital/Orbit.hpp"
#include "gnc/math/Constants.hpp"

#include <cmath>

namespace gnc::orbital {

    Orbit::Orbit(const OrbitalElements& elements)
        : elements_(elements),
          epoch_time_(elements.epoch_time)
    {
    }

    double Orbit::semiMajorAxis() const
    {
        return elements_.semi_major_axis;
    }

    double Orbit::perigeeRadius() const
    {
        return elements_.semi_major_axis *
               (1.0 - elements_.eccentricity);
    }

    double Orbit::apogeeRadius() const
    {
        return elements_.semi_major_axis *
               (1.0 + elements_.eccentricity);
    }

    double Orbit::apogeeAltitude() const
    {
        return apogeeRadius() - gnc::constants::r_E;
    }

    double Orbit::period() const
    {
        return 2.0 * gnc::constants::pi *
               std::sqrt(
                   std::pow(elements_.semi_major_axis, 3) /
                   gnc::constants::mu_E
               );
    }

    double Orbit::meanMotion() const
    {
        return std::sqrt(
            gnc::constants::mu_E /
            std::pow(elements_.semi_major_axis, 3)
        );
    }

    double Orbit::meanAnomaly(double time) const
    {
        return meanMotion() * (time - epoch_time_);
    }

    double Orbit::eccentricAnomaly(double time) const
    {
        double M = meanAnomaly(time);
        double e = elements_.eccentricity;

        double E = M;

        for (int i = 0; i < 100; ++i)
        {
            double f =
                E - e * std::sin(E) - M;

            double df =
                1.0 - e * std::cos(E);

            double delta = f / df;

            E -= delta;

            if (std::abs(delta) <
                gnc::constants::num_tolerance)
            {
                break;
            }
        }

        return E;
    }

    double Orbit::trueAnomaly(double time) const
    {
        double E = eccentricAnomaly(time);
        double e = elements_.eccentricity;

        return 2.0 * std::atan2(
            std::sqrt(1.0 + e) *
                std::sin(E / 2.0),
            std::sqrt(1.0 - e) *
                std::cos(E / 2.0)
        );
    }

    double Orbit::radius(double time) const
    {
        double E = eccentricAnomaly(time);

        return elements_.semi_major_axis *
               (1.0 -
                elements_.eccentricity *
                std::cos(E));
    }

    double Orbit::velocity(double time) const
    {
        double r = radius(time);

        return std::sqrt(
            gnc::constants::mu_E *
            (
                2.0 / r -
                1.0 / elements_.semi_major_axis
            )
        );
    }

    const OrbitalElements& Orbit::elements() const
    {
        return elements_;
    }

}