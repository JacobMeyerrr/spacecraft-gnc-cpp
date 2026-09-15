#include "gnc/control/PIDController.hpp"
#include "gnc/control/ControlSystem.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"
#include "gnc/math/Polynomial.hpp"
#include "gnc/models/MomentumBiasedSpacecraft.hpp"
#include "gnc/simulation/Integrator.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using gnc::math::Matrix4;
using gnc::math::Vector2;
using gnc::math::Vector4;

constexpr double PI = gnc::constants::pi;
constexpr double DEG2RAD = gnc::constants::d2r;
constexpr double RAD2DEG = gnc::constants::r2d;

constexpr double Ix = 2500.0;
constexpr double Iy = 470.0;
constexpr double Iz = 2500.0;

constexpr double phi_max = 0.06 * DEG2RAD;
constexpr double theta_max = 0.06 * DEG2RAD;
constexpr double psi_max = 0.50 * DEG2RAD;

constexpr double To = 86400.0;
constexpr double w0 = 2.0 * PI / To;

constexpr double TxDC = 3.0e-5;
constexpr double TxSinAmp = 6.0e-5;
constexpr double TxPeakAbs = 9.0e-5;
constexpr double TyAmp = 1.2e-4;
constexpr double TzAmp = 6.0e-4;

const gnc::models::MomentumBiasedSpacecraft SATELLITE(
    Ix,
    Iy,
    Iz,
    To);

struct Controller
{
    std::string id;
    double HB{0.0};
    double K{0.0};
    double tau{0.0};
    double k{0.0};
    std::vector<double> K_trials;
    double Tx_step{0.0};
    double Tx_design{0.0};

    gnc::control::MomentumBiasedDesign parameters() const
    {
        gnc::control::MomentumBiasedDesign d;

        d.pitch_kp =
            2.0 * TyAmp / theta_max;

        d.pitch_tau =
            2.0 / std::sqrt(d.pitch_kp / Iy);

        d.pitch_kd =
            d.pitch_kp * d.pitch_tau;

        d.momentum_bias = HB;
        d.roll_gain = K;
        d.roll_tau = tau;
        d.yaw_gain = k;

        return d;
    }
};

struct Series2
{
    std::vector<double> t;
    std::vector<Vector2> y;
    std::vector<Vector2> rate;
};

struct Series1
{
    std::vector<double> t;
    std::vector<double> y;
    std::vector<double> rate;
};

struct TransferSeries
{
    std::vector<double> t;
    std::vector<double> reference;
    std::vector<double> response;
};

Controller makeController(const std::string& variant)
{
    Controller c;

    if(variant == "design_b")
    {
        c.id = "design_b";
        c.HB = 1200.0;
        c.K = 20.0;
        c.tau = 2.0 / std::sqrt(c.K / Ix);
        c.k = 2.0 * std::sqrt(w0 * Iz / c.HB);
        c.K_trials = {
            2.0 * TxPeakAbs / phi_max,
            2.0,
            10.0,
            20.0
        };
        c.Tx_step = TxSinAmp;
        c.Tx_design = TxPeakAbs;
        return c;
    }

    c.id = "design_a";
    c.HB = 2000.0;
    c.K = 12.5;
    c.tau = 2.0 / std::sqrt(c.K / Ix);
    c.k = 2.0 * std::sqrt(w0 * Iz / c.HB);
    c.K_trials = {
        2.0 * TxDC / phi_max,
        1.0,
        5.0,
        12.5
    };
    c.Tx_step = TxDC;
    c.Tx_design = TxDC;
    return c;
}

Matrix4 rollYawOpenA(const Controller& c)
{
    return SATELLITE.openLoopRollYaw(c.HB).A;
}

Eigen::Matrix<double, 4, 2> rollYawB(const Controller& c)
{
    return SATELLITE.openLoopRollYaw(c.HB).B;
}

Matrix4 rollYawClosedA(const Controller& c)
{
    return SATELLITE.closedLoopRollYaw(c.parameters()).A;
}

// ---------------------------------------------------------------
// Roll/yaw constant-input simulation.
// ---------------------------------------------------------------

Series2 simulateRollYawStep(
    const Controller &c,
    bool closed_loop,
    double t_final,
    double dt,
    const Vector2 &input)
{
    Series2 data;
    Vector4 x = Vector4::Zero();

    const Matrix4 A = closed_loop
        ? rollYawClosedA(c)
        : rollYawOpenA(c);

    const auto B = rollYawB(c);

    const auto rhs = [&](double, const Vector4 &state)
    {
        return A * state + B * input;
    };

    const int N =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    data.t.reserve(N);
    data.y.reserve(N);
    data.rate.reserve(N);

    for(int k = 0; k < N; ++k)
    {
        const double t =
            std::min(t_final, static_cast<double>(k) * dt);

        data.t.push_back(t);
        data.y.push_back(
            Vector2(x(0), x(2))
        );
        data.rate.push_back(
            Vector2(x(1), x(3))
        );

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        x = gnc::simulation::rk4Step<Vector4>(
            rhs,
            t,
            x,
            h
        );
    }

    return data;
}

// ---------------------------------------------------------------
// Pitch constant-input simulation.
// ---------------------------------------------------------------

Series1 simulatePitchStep(
    bool closed_loop,
    double t_final,
    double dt,
    double disturbance)
{
    // Pitch design is identical for both design_a and design_b.
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double taup =
        2.0 / wn;

    const double Kd =
        Kp * taup;

    Series1 data;
    Eigen::Vector2d x =
        Eigen::Vector2d::Zero();

    const int N =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    data.t.reserve(N);
    data.y.reserve(N);
    data.rate.reserve(N);

    const auto rhs = [&](double, const Eigen::Vector2d &state)
    {
        const double theta = state(0);
        const double theta_dot = state(1);

        const double control =
            closed_loop
            ? -Kp * theta - Kd * theta_dot
            : 0.0;

        return Eigen::Vector2d(
            theta_dot,
            (control + disturbance) / Iy
        );
    };

    for(int k = 0; k < N; ++k)
    {
        const double t =
            std::min(t_final, static_cast<double>(k) * dt);

        data.t.push_back(t);
        data.y.push_back(x(0));
        data.rate.push_back(x(1));

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        x = gnc::simulation::rk4Step<Eigen::Vector2d>(
            rhs,
            t,
            x,
            h
        );
    }

    return data;
}

// ---------------------------------------------------------------
// Periodic disturbances exactly as specified.
// ---------------------------------------------------------------

Vector2 periodicRollYawDisturbance(double t)
{
    return Vector2(
        TxDC * (1.0 - 2.0 * std::sin(w0 * t)),
        -TzAmp * std::cos(w0 * t)
    );
}

double periodicPitchDisturbance(double t)
{
    return TyAmp * std::cos(w0 * t);
}

struct RollYawPeriodicState
{
    double phi{0.0};
    double phi_dot{0.0};
    double psi{0.0};
    double psi_dot{0.0};
};

RollYawPeriodicState rollYawPeriodicDerivative(
    const Controller &c,
    bool closed_loop,
    double t,
    const RollYawPeriodicState &state)
{
    const double HB = c.HB;

    double roll_stiffness = 0.0;
    double roll_damping = 0.0;
    double yaw_phi_coupling = 0.0;
    double yaw_phidot_coupling = 0.0;

    if(closed_loop)
    {
        roll_stiffness =
            -(w0 * HB + c.K) / Ix;

        roll_damping =
            -(c.K * c.tau) / Ix;

        yaw_phi_coupling =
            c.k * c.K / Iz;

        yaw_phidot_coupling =
            c.k * c.K * c.tau / Iz;
    }
    else
    {
        roll_stiffness =
            -(w0 * HB) / Ix;

        roll_damping = 0.0;

        yaw_phi_coupling = 0.0;
        yaw_phidot_coupling =
            HB / Iz;
    }

    const double disturbance_x =
        TxDC * (1.0 - 2.0 * std::sin(w0 * t));

    const double disturbance_z =
        -TzAmp * std::cos(w0 * t);

    RollYawPeriodicState derivative;

    // These equations are exactly the Lecture-18 linearized
    // momentum-biased roll-yaw equations used in the MATLAB
    // verification: x = [phi, phi_dot, psi, psi_dot]^T.
    derivative.phi =
        state.phi_dot;

    derivative.phi_dot =
        roll_stiffness * state.phi
        + roll_damping * state.phi_dot
        - (HB / Ix) * state.psi_dot
        + disturbance_x / Ix;

    derivative.psi =
        state.psi_dot;

    derivative.psi_dot =
        yaw_phi_coupling * state.phi
        + yaw_phidot_coupling * state.phi_dot
        - (w0 * HB / Iz) * state.psi
        + disturbance_z / Iz;

    return derivative;
}

RollYawPeriodicState addRollYawPeriodicState(
    const RollYawPeriodicState &a,
    const RollYawPeriodicState &b,
    double scale)
{
    RollYawPeriodicState result;

    result.phi =
        a.phi + scale * b.phi;

    result.phi_dot =
        a.phi_dot + scale * b.phi_dot;

    result.psi =
        a.psi + scale * b.psi;

    result.psi_dot =
        a.psi_dot + scale * b.psi_dot;

    return result;
}

Series2 simulateRollYawPeriodic(
    const Controller &c,
    bool closed_loop,
    double t_final,
    double dt)
{
    Series2 data;

    // Use four explicit scalar state variables here instead of
    // storing the periodically forced trajectory in a vector of
    // Eigen fixed-size objects. This makes the validation path
    // deterministic and keeps the implementation visibly aligned
    // with the Lecture-18 state equations.
    RollYawPeriodicState x;

    const int N =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    data.t.reserve(N);
    data.y.reserve(N);
    data.rate.reserve(N);

    for(int step = 0; step < N; ++step)
    {
        const double t =
            std::min(
                t_final,
                static_cast<double>(step) * dt
            );

        data.t.push_back(t);
        data.y.push_back(
            Vector2(x.phi, x.psi)
        );
        data.rate.push_back(
            Vector2(x.phi_dot, x.psi_dot)
        );

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        // Explicit RK4 for the scalar state.
        const RollYawPeriodicState k1 =
            rollYawPeriodicDerivative(
                c,
                closed_loop,
                t,
                x
            );

        const RollYawPeriodicState k2 =
            rollYawPeriodicDerivative(
                c,
                closed_loop,
                t + 0.5 * h,
                addRollYawPeriodicState(x, k1, 0.5 * h)
            );

        const RollYawPeriodicState k3 =
            rollYawPeriodicDerivative(
                c,
                closed_loop,
                t + 0.5 * h,
                addRollYawPeriodicState(x, k2, 0.5 * h)
            );

        const RollYawPeriodicState k4 =
            rollYawPeriodicDerivative(
                c,
                closed_loop,
                t + h,
                addRollYawPeriodicState(x, k3, h)
            );

        RollYawPeriodicState increment;

        increment.phi =
            (h / 6.0)
            * (
                k1.phi
                + 2.0 * k2.phi
                + 2.0 * k3.phi
                + k4.phi
            );

        increment.phi_dot =
            (h / 6.0)
            * (
                k1.phi_dot
                + 2.0 * k2.phi_dot
                + 2.0 * k3.phi_dot
                + k4.phi_dot
            );

        increment.psi =
            (h / 6.0)
            * (
                k1.psi
                + 2.0 * k2.psi
                + 2.0 * k3.psi
                + k4.psi
            );

        increment.psi_dot =
            (h / 6.0)
            * (
                k1.psi_dot
                + 2.0 * k2.psi_dot
                + 2.0 * k3.psi_dot
                + k4.psi_dot
            );

        x = addRollYawPeriodicState(
            x,
            increment,
            1.0
        );
    }

    return data;
}

Series1 simulatePitchPeriodic(
    bool closed_loop,
    double t_final,
    double dt)
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double taup =
        2.0 / wn;

    const double Kd =
        Kp * taup;

    Series1 data;
    Eigen::Vector2d x =
        Eigen::Vector2d::Zero();

    const int N =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    data.t.reserve(N);
    data.y.reserve(N);
    data.rate.reserve(N);

    const auto rhs = [&](double t, const Eigen::Vector2d &state)
    {
        const double theta = state(0);
        const double theta_dot = state(1);

        const double control =
            closed_loop
            ? -Kp * theta - Kd * theta_dot
            : 0.0;

        return Eigen::Vector2d(
            theta_dot,
            (control + periodicPitchDisturbance(t)) / Iy
        );
    };

    for(int k = 0; k < N; ++k)
    {
        const double t =
            std::min(t_final, static_cast<double>(k) * dt);

        data.t.push_back(t);
        data.y.push_back(x(0));
        data.rate.push_back(x(1));

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        x = gnc::simulation::rk4Step<Eigen::Vector2d>(
            rhs,
            t,
            x,
            h
        );
    }

    return data;
}

// ---------------------------------------------------------------
// Controller effort calculations.
// ---------------------------------------------------------------

Vector2 rollYawEffort(
    const Controller& c,
    double phi,
    double phi_dot,
    double phi_ref = 0.0,
    double phi_ref_dot = 0.0)
{
    const auto effort =
        gnc::control::rollYawControlEffort(
            c.parameters(),
            phi,
            phi_dot,
            phi_ref,
            phi_ref_dot);

    return Vector2(
        effort.first,
        effort.second
    );
}

std::vector<double> rollYawDelta(
    const Controller& c,
    double K_value,
    double tau_value)
{
    auto parameters = c.parameters();
    parameters.roll_gain = K_value;
    parameters.roll_tau = tau_value;
    return SATELLITE.characteristicPolynomial(parameters);
}

// ---------------------------------------------------------------
// Generic transfer-function step simulation.
//
// Given:
//   numerator   = [b_n ... b_0]
//   denominator = [a_n ... a_0]
//
// we construct a controllable companion realization and apply a
// constant step input. This is useful for the HW4 angle-command
// checks because the PD controller contains a zero, so the command
// response includes the derivative-of-reference term.
// ---------------------------------------------------------------

TransferSeries simulateTransferStep(
    std::vector<double> numerator,
    std::vector<double> denominator,
    double step_amplitude,
    double t_final,
    double dt)
{
    const int n =
        static_cast<int>(denominator.size()) - 1;

    const double leading = denominator.front();

    for(double &v : numerator)
    {
        v /= leading;
    }

    for(std::size_t i = 0; i < denominator.size(); ++i)
    {
        denominator[i] /= leading;
    }

    // Pad the numerator with leading zeros so that it has the same
    // length as the denominator.
    if(numerator.size() < denominator.size())
    {
        numerator.insert(
            numerator.begin(),
            denominator.size() - numerator.size(),
            0.0
        );
    }

    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(n, n);

    Eigen::VectorXd B =
        Eigen::VectorXd::Zero(n);

    Eigen::RowVectorXd C =
        Eigen::RowVectorXd::Zero(n);

    for(int i = 1; i < n; ++i)
    {
        A(i - 1, i) = 1.0;
    }

    // Last row: [-a0, -a1, ..., -a_(n-1)]
    for(int j = 0; j < n; ++j)
    {
        const double a_j =
            denominator[n - j];

        A(n - 1, j) =
            -a_j;
    }

    B(n - 1) = 1.0;

    // C = [b0, b1, ..., b_(n-1)]
    for(int j = 0; j < n; ++j)
    {
        C(j) =
            numerator[n - j];
    }

    Eigen::VectorXd x =
        Eigen::VectorXd::Zero(n);

    TransferSeries result;

    const int N =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    result.t.reserve(N);
    result.reference.reserve(N);
    result.response.reserve(N);

    const auto rhs = [&](const Eigen::VectorXd &state)
    {
        return A * state + B * step_amplitude;
    };

    for(int k = 0; k < N; ++k)
    {
        const double t =
            std::min(t_final, static_cast<double>(k) * dt);

        result.t.push_back(t);
        result.reference.push_back(step_amplitude);
        result.response.push_back(
            C.dot(x)
        );

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        const Eigen::VectorXd k1 = rhs(x);
        const Eigen::VectorXd k2 =
            rhs(x + 0.5 * h * k1);
        const Eigen::VectorXd k3 =
            rhs(x + 0.5 * h * k2);
        const Eigen::VectorXd k4 =
            rhs(x + h * k3);

        x +=
            (h / 6.0)
            * (k1 + 2.0*k2 + 2.0*k3 + k4);
    }

    return result;
}

// ---------------------------------------------------------------
// CSV writers
// ---------------------------------------------------------------

void writeRollYawStepCsv(
    const std::string &filename,
    const Series2 &open_data,
    const Series2 &closed_data,
    const Controller &c)
{
    std::ofstream out(filename);
    out << std::setprecision(15);

    out <<
        "time_s,phi_open_deg,psi_open_deg,"
        "phi_closed_deg,psi_closed_deg,Tcx_Nm,Tcz_Nm\n";

    for(std::size_t i = 0; i < open_data.t.size(); ++i)
    {
        const Vector2 closed_effort =
            rollYawEffort(
                c,
                closed_data.y[i](0),
                closed_data.rate[i](0)
            );

        out
            << open_data.t[i] << ","
            << open_data.y[i](0) * RAD2DEG << ","
            << open_data.y[i](1) * RAD2DEG << ","
            << closed_data.y[i](0) * RAD2DEG << ","
            << closed_data.y[i](1) * RAD2DEG << ","
            << closed_effort(0) << ","
            << closed_effort(1) << "\n";
    }
}

void writePitchStepCsv(
    const std::string &filename,
    const Series1 &open_data,
    const Series1 &closed_data)
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double Kd =
        Kp * (2.0 / wn);

    std::ofstream out(filename);
    out << std::setprecision(15);

    out <<
        "time_s,theta_open_deg,theta_closed_deg,Tcy_Nm\n";

    for(std::size_t i = 0; i < open_data.t.size(); ++i)
    {
        const double effort =
            -Kp * closed_data.y[i]
            -Kd * closed_data.rate[i];

        out
            << open_data.t[i] << ","
            << open_data.y[i] * RAD2DEG << ","
            << closed_data.y[i] * RAD2DEG << ","
            << effort << "\n";
    }
}

void writePeriodicCsv(
    const std::string &filename,
    const Series2 &ry_open,
    const Series2 &ry_closed,
    const Series1 &pitch_open,
    const Series1 &pitch_closed,
    const Controller &c)
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double Kd =
        Kp * (2.0 / wn);

    std::ofstream out(filename);
    out << std::setprecision(15);

    out <<
        "time_s,phi_open_deg,phi_closed_deg,"
        "theta_open_deg,theta_closed_deg,"
        "psi_open_deg,psi_closed_deg,"
        "Tcx_Nm,Tcy_Nm,Tcz_Nm\n";

    for(std::size_t i = 0; i < ry_open.t.size(); ++i)
    {
        const Vector2 effortRY =
            rollYawEffort(
                c,
                ry_closed.y[i](0),
                ry_closed.rate[i](0)
            );

        const double effortP =
            -Kp * pitch_closed.y[i]
            -Kd * pitch_closed.rate[i];

        out
            << ry_open.t[i] << ","
            << ry_open.y[i](0) * RAD2DEG << ","
            << ry_closed.y[i](0) * RAD2DEG << ","
            << pitch_open.y[i] * RAD2DEG << ","
            << pitch_closed.y[i] * RAD2DEG << ","
            << ry_open.y[i](1) * RAD2DEG << ","
            << ry_closed.y[i](1) * RAD2DEG << ","
            << effortRY(0) << ","
            << effortP << ","
            << effortRY(1) << "\n";
    }
}

void writeCommandCsv(
    const std::string &filename,
    const TransferSeries &command,
    const std::vector<double> &effort1,
    const std::vector<double> &effort2 = {})
{
    std::ofstream out(filename);
    out << std::setprecision(15);

    out << "time_s,reference_deg,response_deg,effort1_Nm,effort2_Nm\n";

    for(std::size_t i = 0; i < command.t.size(); ++i)
    {
        const double effort2_value =
            effort2.empty()
            ? 0.0
            : effort2[i];

        out
            << command.t[i] << ","
            << command.reference[i] * RAD2DEG << ","
            << command.response[i] * RAD2DEG << ","
            << effort1[i] << ","
            << effort2_value << "\n";
    }
}

// ---------------------------------------------------------------
// Utility for maxima.
// ---------------------------------------------------------------

template <typename Container>
double maxAbs(const Container &values)
{
    double result = 0.0;

    for(const double value : values)
    {
        result =
            std::max(result, std::abs(value));
    }

    return result;
}

// ---------------------------------------------------------------
// Root-locus CSVs
// ---------------------------------------------------------------

void writePitchRootLocus(
    const Controller &c)
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double taup =
        2.0 / wn;

    std::ofstream out(
        "data/validation/hw4/hw4_"
        + c.id
        + "_pitch_root_locus.csv"
    );

    out << "K,real1,imag1,real2,imag2\n";

    const int N = 700;
    const double K_min = 1.0e-4;
    const double K_max = 1.0;

    for(int i = 0; i < N; ++i)
    {
        const double alpha =
            static_cast<double>(i) / (N - 1);

        const double K =
            K_min * std::pow(K_max / K_min, alpha);

        const std::complex<double> discriminant(
            K*K*taup*taup - 4.0 * Iy * K,
            0.0
        );

        const std::complex<double> root1 =
            (-K * taup + std::sqrt(discriminant))
            / (2.0 * Iy);

        const std::complex<double> root2 =
            (-K * taup - std::sqrt(discriminant))
            / (2.0 * Iy);

        out
            << K << ","
            << root1.real() << "," << root1.imag() << ","
            << root2.real() << "," << root2.imag() << "\n";
    }
}

void writeRollYawRootLocus(
    const Controller &c)
{
    std::ofstream out(
        "data/validation/hw4/hw4_"
        + c.id
        + "_rollyaw_root_locus.csv"
    );

    out <<
        "K,r1_real,r1_imag,r2_real,r2_imag,"
        "r3_real,r3_imag,r4_real,r4_imag\n";

    const int N = 700;
    const double K_min = 1.0e-4;
    const double K_max = 40.0;

    for(int i = 0; i < N; ++i)
    {
        const double alpha =
            static_cast<double>(i) / (N - 1);

        const double K =
            K_min * std::pow(K_max / K_min, alpha);

        // MATLAB source uses the FINAL tau for the root-locus
        // transfer function, while the design iteration table
        // recomputes tau for each candidate K.
        const auto coeff =
            rollYawDelta(
                c,
                K,
                c.tau
            );

        const auto roots =
            gnc::math::polynomialRoots(coeff);

        out << K;

        for(const auto &root : roots)
        {
            out
                << "," << root.real()
                << "," << root.imag();
        }

        out << "\n";
    }
}

// ---------------------------------------------------------------
// Pitch Bode data and margin calculation
// ---------------------------------------------------------------

void writePitchBode(
    const Controller &c)
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double taup =
        2.0 / wn;

    std::ofstream out(
        "data/validation/hw4/hw4_"
        + c.id
        + "_pitch_bode.csv"
    );

    out << std::setprecision(15);
    out << "omega_rad_s,magnitude,phase_deg\n";

    const int N = 1200;
    const double w_min = 1.0e-5;
    const double w_max = 10.0;

    for(int i = 0; i < N; ++i)
    {
        const double alpha =
            static_cast<double>(i) / (N - 1);

        const double w =
            w_min * std::pow(w_max / w_min, alpha);

        const std::complex<double> s(0.0, w);

        const std::complex<double> L =
            Kp * (taup * s + 1.0)
            / (Iy * s * s);

        double phase =
            std::arg(L) * RAD2DEG;

        // Keep the phase on the conventional -180 to -90 degree
        // branch for this double-integrator plus zero loop.
        if(phase > 0.0)
        {
            phase -= 360.0;
        }

        out
            << w << ","
            << std::abs(L) << ","
            << phase << "\n";
    }
}

double pitchGainCrossover()
{
    const double Kp =
        2.0 * TyAmp / theta_max;

    const double wn =
        std::sqrt(Kp / Iy);

    const double taup =
        2.0 / wn;

    const auto magnitude = [&](double w)
    {
        const std::complex<double> s(0.0, w);

        return std::abs(
            Kp * (taup * s + 1.0)
            / (Iy * s * s)
        );
    };

    double lo = 1.0e-8;
    double hi = 1.0;

    for(int i = 0; i < 100; ++i)
    {
        const double mid =
            0.5 * (lo + hi);

        if(magnitude(mid) > 1.0)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return 0.5 * (lo + hi);
}

// ---------------------------------------------------------------
// Print controller design + verification.
// ---------------------------------------------------------------

void runController(const Controller &c)
{
    // ------------------------
    // Pitch design
    // ------------------------

    const double KpMin =
        TyAmp / theta_max;

    const gnc::control::PDGains pitch =
        gnc::control::designSecondOrderPD(
            Iy,
            TyAmp,
            theta_max,
            2.0,
            1.0);

    const double Kp = pitch.kp;
    const double wnP = pitch.wn;
    const double tauP = pitch.tau;
    const double Kd = pitch.kd;

    // ------------------------
    // Roll/yaw design
    // ------------------------

    const double HB_estimate =
        TzAmp / (w0 * psi_max);

    const double HB_exact =
        (TzAmp + c.k * c.Tx_design)
        / (w0 * psi_max);

    const double HB_large_min =
        std::max({
            Ix * w0,
            Iy * w0,
            Iz * w0
        });

    const double Kmin =
        c.Tx_design / phi_max;

    const double factorization_ratio =
        c.K * c.tau * Iz
        / (c.k * c.HB * Ix);

    std::cout
        << "\n============================================================\n"
        << "HW4 CONTROLLER: " << c.id << "\n"
        << "============================================================\n\n";

    std::cout
        << "Pitch design\n"
        << "Kp_min = " << KpMin << " N*m/rad\n"
        << "Kp     = " << Kp << " N*m/rad\n"
        << "wn     = " << wnP << " rad/s\n"
        << "tau_p  = " << tauP << " s\n"
        << "Kd     = " << Kd << " N*m*s/rad\n"
        << "theta_ss = "
        << (TyAmp / Kp) * RAD2DEG
        << " deg\n\n";

    std::cout
        << "Roll-yaw design\n"
        << "HB estimate = " << HB_estimate << " N*m*s\n"
        << "HB selected = " << c.HB << " N*m*s\n"
        << "HB exact required = " << HB_exact << " N*m*s\n"
        << "HB / max(I*w0) = "
        << c.HB / HB_large_min << "\n"
        << "K_min = " << Kmin << " N*m/rad\n"
        << "K = " << c.K << " N*m/rad\n"
        << "tau = " << c.tau << " s\n"
        << "k = " << c.k << "\n"
        << "factorization ratio = "
        << factorization_ratio << "\n\n";

    // ------------------------
    // Design iterations
    // ------------------------

    std::cout
        << "Roll-yaw K design iterations:\n"
        << "K\t tau [s]\t ratio\n";

    for(const double Ktrial : c.K_trials)
    {
        const double tau_trial =
            2.0 / std::sqrt(Ktrial / Ix);

        const double ratio_trial =
            Ktrial * tau_trial * Iz
            / (c.k * c.HB * Ix);

        std::cout
            << Ktrial << "\t"
            << tau_trial << "\t"
            << ratio_trial << "\n";
    }

    // ------------------------
    // Exact poles of the state-space model
    // ------------------------

    const Matrix4 Acl =
        rollYawClosedA(c);

    Eigen::EigenSolver<Matrix4> pole_solver(Acl);

    std::cout
        << "\nExact closed-loop roll-yaw state-space poles:\n";

    for(int i = 0; i < 4; ++i)
    {
        const auto p =
            pole_solver.eigenvalues()(i);

        std::cout
            << "  "
            << p.real()
            << " "
            << (p.imag() >= 0.0 ? "+ " : "- ")
            << std::abs(p.imag())
            << "j\n";
    }

    const auto delta_roots =
        gnc::math::polynomialRoots(
            rollYawDelta(c, c.K, c.tau)
        );

    std::cout
        << "\nRoots of the Lecture-18 Delta(s) polynomial:\n";

    for(const auto &p : delta_roots)
    {
        std::cout
            << "  "
            << p.real()
            << " "
            << (p.imag() >= 0.0 ? "+ " : "- ")
            << std::abs(p.imag())
            << "j\n";
    }

    // ------------------------
    // Root-locus / Bode CSV data
    // ------------------------

    writePitchRootLocus(c);
    writeRollYawRootLocus(c);
    writePitchBode(c);

    // ------------------------
    // Problem 2: step disturbance
    // ------------------------

    const double ry_step_time =
        To;

    const double pitch_step_time =
        5000.0;

    const double step_dt =
        1.0;

    const Series2 x_open =
        simulateRollYawStep(
            c,
            false,
            ry_step_time,
            step_dt,
            Vector2(c.Tx_step, 0.0)
        );

    const Series2 x_closed =
        simulateRollYawStep(
            c,
            true,
            ry_step_time,
            step_dt,
            Vector2(c.Tx_step, 0.0)
        );

    const Series2 z_open =
        simulateRollYawStep(
            c,
            false,
            ry_step_time,
            step_dt,
            Vector2(0.0, -TzAmp)
        );

    const Series2 z_closed =
        simulateRollYawStep(
            c,
            true,
            ry_step_time,
            step_dt,
            Vector2(0.0, -TzAmp)
        );

    const Series1 pitch_open =
        simulatePitchStep(
            false,
            pitch_step_time,
            step_dt,
            TyAmp
        );

    const Series1 pitch_closed =
        simulatePitchStep(
            true,
            pitch_step_time,
            step_dt,
            TyAmp
        );

    writeRollYawStepCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_roll_step.csv",
        x_open,
        x_closed,
        c
    );

    writeRollYawStepCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_yaw_step.csv",
        z_open,
        z_closed,
        c
    );

    writePitchStepCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_pitch_step.csv",
        pitch_open,
        pitch_closed
    );

    auto maxComponent = [](
        const Series2 &data,
        int component)
    {
        double result = 0.0;

        for(const auto &v : data.y)
        {
            result =
                std::max(
                    result,
                    std::abs(v(component))
                );
        }

        return result;
    };

    std::cout
        << "\nProblem 2 step disturbance maxima [deg]:\n"
        << "x-step magnitude = "
        << c.Tx_step
        << " N*m\n"
        << "x-step open phi/psi = "
        << maxComponent(x_open, 0) * RAD2DEG
        << " / "
        << maxComponent(x_open, 1) * RAD2DEG
        << "\n"
        << "x-step closed phi/psi = "
        << maxComponent(x_closed, 0) * RAD2DEG
        << " / "
        << maxComponent(x_closed, 1) * RAD2DEG
        << "\n"
        << "y-step closed pitch = "
        << maxAbs(pitch_closed.y) * RAD2DEG
        << "\n"
        << "z-step open phi/psi = "
        << maxComponent(z_open, 0) * RAD2DEG
        << " / "
        << maxComponent(z_open, 1) * RAD2DEG
        << "\n"
        << "z-step closed phi/psi = "
        << maxComponent(z_closed, 0) * RAD2DEG
        << " / "
        << maxComponent(z_closed, 1) * RAD2DEG
        << "\n";

    // ------------------------
    // Problem 3: periodic disturbance
    // ------------------------

    const double periodic_time =
        2.0 * To;

    const double periodic_dt =
        1.0;

    const Series2 periodic_ry_open =
        simulateRollYawPeriodic(
            c,
            false,
            periodic_time,
            periodic_dt
        );

    const Series2 periodic_ry_closed =
        simulateRollYawPeriodic(
            c,
            true,
            periodic_time,
            periodic_dt
        );

    const Series1 periodic_pitch_open =
        simulatePitchPeriodic(
            false,
            periodic_time,
            periodic_dt
        );

    const Series1 periodic_pitch_closed =
        simulatePitchPeriodic(
            true,
            periodic_time,
            periodic_dt
        );

    writePeriodicCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_periodic.csv",
        periodic_ry_open,
        periodic_ry_closed,
        periodic_pitch_open,
        periodic_pitch_closed,
        c
    );

    const double periodic_open_roll_max =
        maxComponent(periodic_ry_open, 0) * RAD2DEG;

    const double periodic_open_pitch_max =
        maxAbs(periodic_pitch_open.y) * RAD2DEG;

    const double periodic_open_yaw_max =
        maxComponent(periodic_ry_open, 1) * RAD2DEG;

    const double periodic_closed_roll_max =
        maxComponent(periodic_ry_closed, 0) * RAD2DEG;

    const double periodic_closed_pitch_max =
        maxAbs(periodic_pitch_closed.y) * RAD2DEG;

    const double periodic_closed_yaw_max =
        maxComponent(periodic_ry_closed, 1) * RAD2DEG;

    std::cout
        << "\nProblem 3 periodic disturbance maxima [deg]:\n"
        << "open roll/pitch/yaw = "
        << periodic_open_roll_max
        << " / "
        << periodic_open_pitch_max
        << " / "
        << periodic_open_yaw_max
        << "\n"
        << "closed roll/pitch/yaw = "
        << periodic_closed_roll_max
        << " / "
        << periodic_closed_pitch_max
        << " / "
        << periodic_closed_yaw_max
        << "\n"
        << "Reference closed-loop targets (Design B): "
        << "0.05087 / 0.03000 / 0.39436 deg\n";

    // ------------------------
    // Problem 3 control effort
    // ------------------------

    double max_Tcx_periodic = 0.0;
    double max_Tcz_periodic = 0.0;
    double max_Tcy_periodic = 0.0;

    for(std::size_t i = 0;
        i < periodic_ry_closed.t.size();
        ++i)
    {
        const Vector2 effortRY =
            rollYawEffort(
                c,
                periodic_ry_closed.y[i](0),
                periodic_ry_closed.rate[i](0)
            );

        const double KpLocal =
            2.0 * TyAmp / theta_max;

        const double taupLocal =
            2.0 / std::sqrt(KpLocal / Iy);

        const double KdLocal =
            KpLocal * taupLocal;

        const double effortP =
            -KpLocal * periodic_pitch_closed.y[i]
            -KdLocal * periodic_pitch_closed.rate[i];

        max_Tcx_periodic =
            std::max(
                max_Tcx_periodic,
                std::abs(effortRY(0))
            );

        max_Tcz_periodic =
            std::max(
                max_Tcz_periodic,
                std::abs(effortRY(1))
            );

        max_Tcy_periodic =
            std::max(
                max_Tcy_periodic,
                std::abs(effortP)
            );
    }

    std::cout
        << "\nPeriodic maximum control effort [N*m]:\n"
        << "max |Tcx| = " << max_Tcx_periodic << "\n"
        << "max |Tcy| = " << max_Tcy_periodic << "\n"
        << "max |Tcz| = " << max_Tcz_periodic << "\n"
        << "Reference Design-B effort: "
        << "1.794e-2 / 1.362e-4 / 4.416e-4 N*m\n";

    // ------------------------
    // Angle-command verification
    // ------------------------

    // Pitch command transfer function from the submitted MATLAB:
    //
    //   theta/theta_ref = (Kd*s + Kp)
    //                    / (Iy*s^2 + Kd*s + Kp)
    const std::vector<double> pitch_num = {
        Kd,
        Kp
    };

    const std::vector<double> pitch_den = {
        Iy,
        Kd,
        Kp
    };

    const TransferSeries pitch_command =
        simulateTransferStep(
            pitch_num,
            pitch_den,
            theta_max,
            5000.0,
            0.1
        );

    std::vector<double> pitch_command_effort(
        pitch_command.t.size(),
        0.0
    );

    for(std::size_t i = 0;
        i < pitch_command.t.size();
        ++i)
    {
        // Reporting convention from the reference MATLAB:
        // control = Kp*(ref - theta)
        //          + Kd*(ref_dot - theta_dot)
        // with ref_dot = 0 for the constant step.
        // Estimate theta_dot with a central difference below.
        double theta_dot = 0.0;

        if(i == 0 && pitch_command.t.size() > 1)
        {
            theta_dot =
                (
                    pitch_command.response[1]
                    - pitch_command.response[0]
                )
                / (
                    pitch_command.t[1]
                    - pitch_command.t[0]
                );
        }
        else if(i + 1 < pitch_command.t.size())
        {
            theta_dot =
                (
                    pitch_command.response[i+1]
                    - pitch_command.response[i-1]
                )
                / (
                    pitch_command.t[i+1]
                    - pitch_command.t[i-1]
                );
        }
        else if(i > 0)
        {
            theta_dot =
                (
                    pitch_command.response[i]
                    - pitch_command.response[i-1]
                )
                / (
                    pitch_command.t[i]
                    - pitch_command.t[i-1]
                );
        }

        pitch_command_effort[i] =
            Kp * (
                theta_max
                - pitch_command.response[i]
            )
            - Kd * theta_dot;
    }

    writeCommandCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_pitch_command.csv",
        pitch_command,
        pitch_command_effort
    );

    // Roll reference transfer function from the reference design
    // Appendix A / Lecture 18 form:
    //
    //   phi/phi_ref =
    //       K(tau s+1)(Iz s^2 + k HB s + w0 HB)
    //       / Delta(s)
    const std::vector<double> roll_num = {
        c.K * c.tau * Iz,
        c.K * (
            Iz
            + c.k * c.tau * c.HB
        ),
        c.K * (
            w0 * c.HB * c.tau
            + c.k * c.HB
        ),
        c.K * w0 * c.HB
    };

    const std::vector<double> roll_den =
        rollYawDelta(
            c,
            c.K,
            c.tau
        );

    const TransferSeries roll_command =
        simulateTransferStep(
            roll_num,
            roll_den,
            phi_max,
            5000.0,
            0.1
        );

    std::vector<double> roll_command_Tcx(
        roll_command.t.size(),
        0.0
    );

    std::vector<double> roll_command_Tcz(
        roll_command.t.size(),
        0.0
    );

    // Numerical derivative for the effort report.
    std::vector<double> phi_dot(
        roll_command.t.size(),
        0.0
    );

    for(std::size_t i = 0;
        i < roll_command.t.size();
        ++i)
    {
        if(i == 0 && roll_command.t.size() > 1)
        {
            phi_dot[i] =
                (
                    roll_command.response[1]
                    - roll_command.response[0]
                )
                / (
                    roll_command.t[1]
                    - roll_command.t[0]
                );
        }
        else if(i + 1 < roll_command.t.size())
        {
            phi_dot[i] =
                (
                    roll_command.response[i+1]
                    - roll_command.response[i-1]
                )
                / (
                    roll_command.t[i+1]
                    - roll_command.t[i-1]
                );
        }
        else if(i > 0)
        {
            phi_dot[i] =
                (
                    roll_command.response[i]
                    - roll_command.response[i-1]
                )
                / (
                    roll_command.t[i]
                    - roll_command.t[i-1]
                );
        }

        const double q =
            c.tau * phi_dot[i]
            + (
                roll_command.response[i]
                - phi_max
            );

        roll_command_Tcx[i] =
            c.K * q;

        roll_command_Tcz[i] =
            -c.k * roll_command_Tcx[i];
    }

    writeCommandCsv(
        "data/validation/hw4/hw4_"
        + c.id
        + "_roll_command.csv",
        roll_command,
        roll_command_Tcx,
        roll_command_Tcz
    );

    std::cout
        << "\nAngle-command verification:\n"
        << "0.06-deg pitch command peak/final = "
        << maxAbs(pitch_command.response) * RAD2DEG
        << " / "
        << pitch_command.response.back() * RAD2DEG
        << " deg\n"
        << "0.06-deg roll command peak/final = "
        << maxAbs(roll_command.response) * RAD2DEG
        << " / "
        << roll_command.response.back() * RAD2DEG
        << " deg\n";

    // ------------------------
    // Pitch frequency response
    // ------------------------

    const double wgc =
        pitchGainCrossover();

    const double phase_margin =
        std::atan2(
            Kd * wgc,
            Kp
        ) * RAD2DEG;

    std::cout
        << "\nPitch frequency-domain verification:\n"
        << "w_gc = " << wgc << " rad/s\n"
        << "phase margin = " << phase_margin << " deg\n"
        << "gain margin = infinity\n\n";
}

}

int main(int argc, char **argv)
{
    std::filesystem::create_directories("data/validation/hw4");
    std::cout << std::setprecision(10);

    std::string variant =
        argc > 1
        ? argv[1]
        : "design_b";

    if(variant == "all")
    {
        runController(
            makeController("design_a")
        );

        runController(
            makeController("design_b")
        );

        std::cout
            << "\nFinished BOTH reference controller variants.\n";

        return 0;
    }

    if(variant != "design_a" && variant != "design_b")
    {
        std::cerr
            << "Usage: ./build/validate_hw4 [design_a|design_b|all]\n";

        return 1;
    }

    runController(
        makeController(variant)
    );

    std::cout
        << "\nFinished HW4 controller variant: "
        << variant
        << "\n";

    return 0;
}
