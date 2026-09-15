#include "gnc/orbital/OrbitalElements.hpp"
#include "gnc/orbital/Orbit.hpp"
#include "gnc/math/Constants.hpp"

#include <iostream>
#include <iomanip>

int main()
{
    // HW1 solution-case inputs
    double perigee_altitude = 500'000.0; // m
    double eccentricity = 0.9;

    double perigee_radius =
        gnc::constants::r_E + perigee_altitude;

    double semi_major_axis =
        perigee_radius / (1.0 - eccentricity);

    gnc::orbital::OrbitalElements elements{
        semi_major_axis,
        eccentricity,
        30.0 * gnc::constants::d2r,
        45.0 * gnc::constants::d2r,
        -30.0 * gnc::constants::d2r,
        0.0
    };

    gnc::orbital::Orbit orbit(elements);

    std::cout << std::setprecision(12);

    std::cout << "Semi-major axis: "
              << orbit.semiMajorAxis() << " m\n";

    std::cout << "Orbital period: "
              << orbit.period() << " s\n";

    std::cout << "Apogee altitude: "
              << orbit.apogeeAltitude() << " m\n";

    std::cout << "Mean motion: "
              << orbit.meanMotion() << " rad/s\n";

    std::cout << "\nTime [s], True Anomaly [rad], Velocity [m/s]\n";

    int samples = 100;

    for (int i = 0; i <= samples; ++i)
    {
        double t =
            orbit.period() *
            static_cast<double>(i) /
            samples;

        std::cout
            << t << ", "
            << orbit.trueAnomaly(t) << ", "
            << orbit.velocity(t) << '\n';
    }

    return 0;
}