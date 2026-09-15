#pragma once

#include "gnc/control/ControlSystem.hpp"

namespace gnc::models {

    class MomentumBiasedSpacecraft
    {
        public:
            MomentumBiasedSpacecraft(
                double ix,
                double iy,
                double iz,
                double orbital_period);

            double Ix() const;
            double Iy() const;
            double Iz() const;
            double orbitalRate() const;
            double orbitalPeriod() const;

            gnc::control::StateSpaceModel openLoopRollYaw(
                double momentum_bias) const;

            gnc::control::StateSpaceModel closedLoopRollYaw(
                const gnc::control::MomentumBiasedDesign& design) const;

            std::vector<double> characteristicPolynomial(
                const gnc::control::MomentumBiasedDesign& design) const;

        private:
            double ix_;
            double iy_;
            double iz_;
            double orbital_period_;
            double omega0_;
    };

}
