#include "gnc/control/PIDController.hpp"

#include <cmath>
#include <stdexcept>

namespace gnc::control {

    PIDController::PIDController(PIDGains gains)
        : gains_(gains)
    {
    }

    void PIDController::setGains(PIDGains gains)
    {
        gains_ = gains;
    }

    void PIDController::reset()
    {
        integral_error_ = 0.0;
    }

    double PIDController::update(
        double error,
        double error_rate,
        double dt)
    {
        if(dt < 0.0)
        {
            throw std::invalid_argument("Controller dt cannot be negative");
        }

        integral_error_ += error * dt;

        return gains_.kp * error
             + gains_.ki * integral_error_
             + gains_.kd * error_rate;
    }

    const PIDGains& PIDController::gains() const
    {
        return gains_;
    }

    PDGains designSecondOrderPD(
        double inertia,
        double disturbance_torque,
        double max_attitude_error_rad,
        double gain_margin_factor,
        double damping_ratio)
    {
        if(inertia <= 0.0 ||
           max_attitude_error_rad <= 0.0 ||
           gain_margin_factor <= 0.0 ||
           damping_ratio < 0.0)
        {
            throw std::invalid_argument(
                "Invalid second-order PD design inputs"
            );
        }

        PDGains gains;

        const double kp_min =
            disturbance_torque
            / max_attitude_error_rad;

        gains.kp =
            gain_margin_factor * kp_min;

        gains.wn =
            std::sqrt(gains.kp / inertia);

        gains.zeta =
            damping_ratio;

        gains.tau =
            2.0 * gains.zeta / gains.wn;

        gains.kd =
            gains.kp * gains.tau;

        return gains;
    }

}
