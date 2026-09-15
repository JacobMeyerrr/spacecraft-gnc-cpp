#pragma once

#include "gnc/orbital/OrbitalElements.hpp"

namespace gnc::orbital {

    class Orbit
    {
        public:
            explicit Orbit(const OrbitalElements& elements);

            double semiMajorAxis() const;
            double perigeeRadius() const;
            double apogeeRadius() const;
            double apogeeAltitude() const;
            double period() const;
            double meanMotion() const;

            // Analytical two-body helpers for HW1
            double meanAnomaly(double time) const;
            double eccentricAnomaly(double time) const;
            double trueAnomaly(double time) const;
            double radius(double time) const;
            double velocity(double time) const;

            const OrbitalElements& elements() const;

        private:
            OrbitalElements elements_;
            double epoch_time_;
    };

}