#pragma once

#include "gnc/control/StateSpace.hpp"

#include <utility>
#include <vector>

namespace gnc::control {

    enum class MomentumBiasedDesignVariant
    {
        DesignA,
        DesignB
    };

    struct MomentumBiasedDesign
    {
        double pitch_kp{0.0};
        double pitch_tau{0.0};
        double pitch_kd{0.0};
        double momentum_bias{0.0};
        double roll_gain{0.0};
        double roll_tau{0.0};
        double yaw_gain{0.0};
    };

    MomentumBiasedDesign designMomentumBiasedController(
        MomentumBiasedDesignVariant variant,
        double ix,
        double iy,
        double iz,
        double orbital_period,
        double roll_limit_rad,
        double pitch_limit_rad,
        double yaw_limit_rad,
        double tx_amplitude,
        double ty_amplitude,
        double tz_amplitude);

    std::vector<double> rollYawCharacteristicPolynomial(
        const MomentumBiasedDesign& design,
        double ix,
        double iz,
        double orbital_rate);

    StateSpaceModel makeRollYawOpenLoop(
        double ix,
        double iz,
        double orbital_rate,
        double momentum_bias);

    StateSpaceModel makeRollYawClosedLoop(
        double ix,
        double iz,
        double orbital_rate,
        const MomentumBiasedDesign& design);

    double pitchGainCrossover(
        double iy,
        double kp,
        double tau);

    double pitchPhaseMargin(
        double iy,
        double kp,
        double tau);

    std::pair<double, double> rollYawControlEffort(
        const MomentumBiasedDesign& design,
        double phi,
        double phi_dot,
        double phi_reference = 0.0,
        double phi_reference_dot = 0.0);

}
