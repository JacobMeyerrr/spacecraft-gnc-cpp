#pragma once

namespace gnc::control {

    struct PDGains
    {
        double kp{0.0};
        double kd{0.0};
        double tau{0.0};
        double wn{0.0};
        double zeta{0.0};
    };

    struct PIDGains
    {
        double kp{0.0};
        double ki{0.0};
        double kd{0.0};
    };

    class PIDController
    {
        public:
            explicit PIDController(PIDGains gains = {});
            void setGains(PIDGains gains);
            void reset();
            double update(double error, double error_rate, double dt);
            const PIDGains& gains() const;

        private:
            PIDGains gains_;
            double integral_error_{0.0};
    };

    PDGains designSecondOrderPD(
        double inertia,
        double disturbance_torque,
        double max_attitude_error_rad,
        double gain_margin_factor,
        double damping_ratio);

}
