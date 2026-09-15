#include "gnc/control/ControlSystem.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Polynomial.hpp"

#include <cmath>
#include <stdexcept>

namespace gnc::control {

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
        double tz_amplitude)
    {
        if(ix <= 0.0 || iy <= 0.0 || iz <= 0.0 || orbital_period <= 0.0)
        {
            throw std::invalid_argument("Invalid spacecraft/controller parameters");
        }

        (void)yaw_limit_rad;

        const double w0 =
            2.0 * gnc::constants::pi / orbital_period;

        const double pitch_kp_min =
            ty_amplitude / pitch_limit_rad;

        MomentumBiasedDesign d;
        d.pitch_kp = 2.0 * pitch_kp_min;
        d.pitch_tau = 2.0 / std::sqrt(d.pitch_kp / iy);
        d.pitch_kd = d.pitch_kp * d.pitch_tau;

        if(variant == MomentumBiasedDesignVariant::DesignB)
        {
            d.momentum_bias = 1200.0;
            d.roll_gain = 20.0;
        }
        else
        {
            d.momentum_bias = 2000.0;
            d.roll_gain = 12.5;
        }

        d.yaw_gain =
            2.0 * std::sqrt(
                w0 * iz / d.momentum_bias
            );

        d.roll_tau =
            2.0 / std::sqrt(d.roll_gain / ix);

        (void)roll_limit_rad;
        (void)tx_amplitude;
        (void)tz_amplitude;

        return d;
    }

    std::vector<double> rollYawCharacteristicPolynomial(
        const MomentumBiasedDesign& d,
        double ix,
        double iz,
        double orbital_rate)
    {
        return {
            ix * iz,
            d.roll_gain * iz * d.roll_tau,
            d.roll_gain * iz
                + ix * orbital_rate * d.momentum_bias
                + d.yaw_gain * d.roll_gain
                    * d.roll_tau * d.momentum_bias,
            orbital_rate * d.momentum_bias
                * d.roll_gain * d.roll_tau
                + d.yaw_gain * d.roll_gain * d.momentum_bias,
            orbital_rate * d.momentum_bias * d.roll_gain
        };
    }

    StateSpaceModel makeRollYawOpenLoop(
        double ix,
        double iz,
        double orbital_rate,
        double momentum_bias)
    {
        StateSpaceModel model;

        model.A = gnc::math::Matrix4::Zero();
        model.A(0,1) = 1.0;
        model.A(1,0) = -orbital_rate * momentum_bias / ix;
        model.A(1,3) = -momentum_bias / ix;
        model.A(2,3) = 1.0;
        model.A(3,1) = momentum_bias / iz;
        model.A(3,2) = -orbital_rate * momentum_bias / iz;

        model.B = Eigen::Matrix<double,4,2>::Zero();
        model.B(1,0) = 1.0 / ix;
        model.B(3,1) = 1.0 / iz;

        model.C = Eigen::Matrix<double, 2, 4>::Zero();
        model.C(0,0) = 1.0;
        model.C(1,2) = 1.0;
        model.D = Eigen::Matrix2d::Zero();

        return model;
    }

    StateSpaceModel makeRollYawClosedLoop(
        double ix,
        double iz,
        double orbital_rate,
        const MomentumBiasedDesign& d)
    {
        StateSpaceModel model =
            makeRollYawOpenLoop(
                ix,
                iz,
                orbital_rate,
                d.momentum_bias
            );

        model.A(1,0) =
            -(orbital_rate * d.momentum_bias + d.roll_gain) / ix;

        model.A(1,1) =
            -(d.roll_gain * d.roll_tau) / ix;

        model.A(3,0) =
            d.yaw_gain * d.roll_gain / iz;

        model.A(3,1) =
            d.yaw_gain * d.roll_gain * d.roll_tau / iz;

        return model;
    }

    double pitchGainCrossover(
        double iy,
        double kp,
        double tau)
    {
        const auto magnitude = [=](double w)
        {
            const double numerator =
                kp * std::sqrt(
                    1.0 + tau * tau * w * w
                );

            return numerator / (iy * w * w);
        };

        double low = 1.0e-9;
        double high = 10.0;

        for(int i = 0; i < 100; ++i)
        {
            const double mid = 0.5 * (low + high);

            if(magnitude(mid) > 1.0)
            {
                low = mid;
            }
            else
            {
                high = mid;
            }
        }

        return 0.5 * (low + high);
    }

    double pitchPhaseMargin(
        double iy,
        double kp,
        double tau)
    {
        const double w_c =
            pitchGainCrossover(iy, kp, tau);

        return 180.0 * gnc::constants::d2r
             - (
                -gnc::constants::pi
                + std::atan2(kp * tau * w_c, kp)
             );
    }

    std::pair<double,double> rollYawControlEffort(
        const MomentumBiasedDesign& d,
        double phi,
        double phi_dot,
        double phi_reference,
        double phi_reference_dot)
    {
        const double q =
            d.roll_tau * (phi_dot - phi_reference_dot)
            + (phi - phi_reference);

        const double tcx =
            d.roll_gain * q;

        const double tcz =
            -d.yaw_gain * tcx;

        return {tcx, tcz};
    }

}
