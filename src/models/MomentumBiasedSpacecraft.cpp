#include "gnc/models/MomentumBiasedSpacecraft.hpp"
#include "gnc/math/Constants.hpp"

namespace gnc::models {

    MomentumBiasedSpacecraft::MomentumBiasedSpacecraft(
        double ix,
        double iy,
        double iz,
        double orbital_period)
        : ix_(ix),
          iy_(iy),
          iz_(iz),
          orbital_period_(orbital_period),
          omega0_(2.0 * gnc::constants::pi / orbital_period)
    {
    }

    double MomentumBiasedSpacecraft::Ix() const { return ix_; }
    double MomentumBiasedSpacecraft::Iy() const { return iy_; }
    double MomentumBiasedSpacecraft::Iz() const { return iz_; }
    double MomentumBiasedSpacecraft::orbitalRate() const { return omega0_; }
    double MomentumBiasedSpacecraft::orbitalPeriod() const { return orbital_period_; }

    gnc::control::StateSpaceModel
    MomentumBiasedSpacecraft::openLoopRollYaw(
        double momentum_bias) const
    {
        return gnc::control::makeRollYawOpenLoop(
            ix_,
            iz_,
            omega0_,
            momentum_bias
        );
    }

    gnc::control::StateSpaceModel
    MomentumBiasedSpacecraft::closedLoopRollYaw(
        const gnc::control::MomentumBiasedDesign& design) const
    {
        return gnc::control::makeRollYawClosedLoop(
            ix_,
            iz_,
            omega0_,
            design
        );
    }

    std::vector<double>
    MomentumBiasedSpacecraft::characteristicPolynomial(
        const gnc::control::MomentumBiasedDesign& design) const
    {
        return gnc::control::rollYawCharacteristicPolynomial(
            design,
            ix_,
            iz_,
            omega0_
        );
    }

}
