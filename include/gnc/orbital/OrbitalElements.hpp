#pragma once

namespace gnc::orbital {

    // the common foundational orbital elements
    struct OrbitalElements
    {
        double semi_major_axis;
        double eccentricity;
        double inclination;
        double raan;
        double argument_of_periapsis;
        double true_anomaly;
        double epoch_time;
    };

}