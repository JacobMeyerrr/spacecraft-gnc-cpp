#include "gnc/attitude/AttitudeKinematics.hpp"
#include "gnc/attitude/EulerAngles.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"
#include "gnc/models/EarthPointingSpacecraft.hpp"
#include "gnc/simulation/DataLogger.hpp"
#include "gnc/simulation/Integrator.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace {

using gnc::math::Matrix3;
using gnc::math::Vector3;
using gnc::math::Vector6;

constexpr double MU = 4.035040e14;
constexpr double RE = 6.3753845e6;
constexpr double ALTITUDE = 700.0e3;

const Matrix3 INERTIA = (Matrix3() <<
    90.0, 2.0, -1.0,
    2.0, 120.0, 1.0,
    -1.0, 1.0, 40.0
).finished();

Vector3 disturbance(
    double t,
    int case_id)
{
    if(case_id != 3)
    {
        return Vector3::Zero();
    }

    return 0.1
        * Vector3(-1.0, 1.0, -1.0)
        * (1.0 + std::sin(3.0 * t));
}

struct SimulationResult
{
    std::vector<double> time;
    std::vector<Vector6> linear;
    std::vector<Vector6> nonlinear;
};

SimulationResult simulateCase(
    const gnc::models::EarthPointingSpacecraft& spacecraft,
    int case_id,
    double t_final,
    double dt)
{
    SimulationResult result;

    const Vector3 eta0_deg =
        case_id == 1 ? Vector3(5.0, 2.0, 1.0) :
        case_id == 2 ? Vector3(1.0, -2.0, 5.0) :
                       Vector3::Zero();

    const Vector3 eta0 =
        eta0_deg * gnc::constants::d2r;

    Vector6 x_linear = Vector6::Zero();
    x_linear.head<3>() = eta0;

    Vector6 x_nonlinear = Vector6::Zero();
    x_nonlinear.head<3>() = eta0;

    const gnc::attitude::EulerAngles angles{
        eta0(0), eta0(1), eta0(2)
    };

    // Zero initial motion relative to the orbit means the body-rate
    // relative to inertial space is the orbital reference rate
    // expressed in the body frame.
    x_nonlinear.tail<3>() =
        gnc::attitude::dcm123(angles)
        * Vector3(0.0, -spacecraft.orbitalRate(), 0.0);

    const int sample_count =
        static_cast<int>(std::ceil(t_final / dt)) + 1;

    result.time.reserve(sample_count);
    result.linear.reserve(sample_count);
    result.nonlinear.reserve(sample_count);

    double t = 0.0;

    while(t <= t_final + 1.0e-12)
    {
        result.time.push_back(t);
        result.linear.push_back(x_linear);
        result.nonlinear.push_back(x_nonlinear);

        if(t >= t_final)
        {
            break;
        }

        const double h =
            std::min(dt, t_final - t);

        x_linear = gnc::simulation::rk4Step<Vector6>(
            [&](double stage_time, const Vector6& state)
            {
                return spacecraft.linearDerivative(
                    state,
                    disturbance(stage_time, case_id));
            },
            t,
            x_linear,
            h);

        x_nonlinear = gnc::simulation::rk4Step<Vector6>(
            [&](double stage_time, const Vector6& state)
            {
                return spacecraft.nonlinearDerivative(
                    stage_time,
                    state,
                    disturbance(stage_time, case_id));
            },
            t,
            x_nonlinear,
            h);

        t += h;
    }

    return result;
}

void writeCaseCsv(
    const SimulationResult& result,
    int case_id)
{
    const std::string path =
        "data/validation/hw3/hw3_case"
        + std::to_string(case_id)
        + ".csv";

    gnc::simulation::CsvLogger log(path);
    log.writeHeader(
        "time_s,"
        "phi_linear_deg,theta_linear_deg,psi_linear_deg,"
        "phi_nonlinear_deg,theta_nonlinear_deg,psi_nonlinear_deg,"
        "phi_error_deg,theta_error_deg,psi_error_deg,"
        "Tdx_Nm,Tdy_Nm,Tdz_Nm"
    );

    for(std::size_t i = 0; i < result.time.size(); ++i)
    {
        const Vector3 linear =
            result.linear[i].head<3>() * gnc::constants::r2d;

        const Vector3 nonlinear =
            result.nonlinear[i].head<3>() * gnc::constants::r2d;

        const Vector3 error =
            nonlinear - linear;

        const Vector3 td =
            disturbance(result.time[i], case_id);

        std::ostringstream line;
        line << std::setprecision(15)
             << result.time[i] << ","
             << linear(0) << "," << linear(1) << "," << linear(2) << ","
             << nonlinear(0) << "," << nonlinear(1) << "," << nonlinear(2) << ","
             << error(0) << "," << error(1) << "," << error(2) << ","
             << td(0) << "," << td(1) << "," << td(2);

        log.writeLine(line.str());
    }
}

void writeTransferFunctionGrid(
    const gnc::models::EarthPointingSpacecraft& spacecraft)
{
    std::ofstream out(
        "data/validation/hw3/hw3_frequency_response.csv"
    );

    out << std::setprecision(15);

    out <<
        "omega_rad_s,"
        "G11,G12,G13,G21,G22,G23,G31,G32,G33\n";

    const Matrix3 I = spacecraft.inertia();
    const Matrix3 D = spacecraft.linearDampingMatrix();
    const Matrix3 K = spacecraft.linearStiffnessMatrix();

    constexpr int points = 800;
    constexpr double omega_min = 1.0e-6;
    constexpr double omega_max = 2.0;

    for(int i = 0; i < points; ++i)
    {
        const double fraction =
            static_cast<double>(i)
            / static_cast<double>(points - 1);

        const double omega =
            omega_min
            * std::pow(
                omega_max / omega_min,
                fraction
            );

        const std::complex<double> s(0.0, omega);

        const Eigen::Matrix3cd Z =
            I.cast<std::complex<double>>() * s * s
            + D.cast<std::complex<double>>() * s
            + K.cast<std::complex<double>>();

        const Eigen::Matrix3cd G = Z.inverse();

        out << omega;

        for(int row = 0; row < 3; ++row)
        {
            for(int column = 0; column < 3; ++column)
            {
                out << "," << std::abs(G(row, column));
            }
        }

        out << '\n';
    }
}

void printCaseSummary(
    const SimulationResult& result,
    int case_id)
{
    Vector3 max_error = Vector3::Zero();
    Vector3 rms_sum = Vector3::Zero();

    for(std::size_t i = 0; i < result.time.size(); ++i)
    {
        const Vector3 error =
            (
                result.nonlinear[i].head<3>()
                - result.linear[i].head<3>()
            ) * gnc::constants::r2d;

        for(int axis = 0; axis < 3; ++axis)
        {
            max_error(axis) =
                std::max(
                    max_error(axis),
                    std::abs(error(axis))
                );

            rms_sum(axis) +=
                error(axis) * error(axis);
        }
    }

    const Vector3 rms =
        (rms_sum / static_cast<double>(result.time.size()))
            .array()
            .sqrt();

    const Vector3 final_linear =
        result.linear.back().head<3>()
        * gnc::constants::r2d;

    const Vector3 final_nonlinear =
        result.nonlinear.back().head<3>()
        * gnc::constants::r2d;

    std::cout
        << "Case " << case_id << '\n'
        << "  max |nonlinear-linear| [deg] = "
        << max_error.transpose() << '\n'
        << "  RMS error [deg] = "
        << rms.transpose() << '\n'
        << "  final linear [deg] = "
        << final_linear.transpose() << '\n'
        << "  final nonlinear [deg] = "
        << final_nonlinear.transpose() << "\n\n";
}

}

int main()
{
    std::filesystem::create_directories(
        "data/validation/hw3"
    );

    const gnc::models::EarthPointingSpacecraft spacecraft(
        MU,
        RE,
        ALTITUDE,
        INERTIA
    );

    std::cout << std::setprecision(12);

    std::cout
        << "HW3 ADCS validation\n\n"
        << "Orbital radius: "
        << spacecraft.orbitalRadius()
        << " m\n"
        << "Orbital rate: "
        << spacecraft.orbitalRate()
        << " rad/s\n"
        << "Orbital period: "
        << spacecraft.orbitalPeriod()
        << " s\n\n";

    std::cout
        << "Linear damping matrix D:\n"
        << spacecraft.linearDampingMatrix()
        << "\n\n";

    std::cout
        << "Linear stiffness matrix K:\n"
        << spacecraft.linearStiffnessMatrix()
        << "\n\n";

    std::cout
        << "Nominal trim torque [N m]:\n"
        << spacecraft.trimTorque()
        << "\n\n";

    writeTransferFunctionGrid(spacecraft);

    for(int case_id = 1; case_id <= 3; ++case_id)
    {
        const double t_final =
            case_id == 3
            ? 20.0
            : 2.0 * spacecraft.orbitalPeriod();

        const double dt =
            case_id == 3
            ? 0.002
            : 0.25;

        std::cout
            << "Running Case " << case_id
            << " (dt = " << dt
            << " s, tf = " << t_final
            << " s)\n";

        const SimulationResult result =
            simulateCase(
                spacecraft,
                case_id,
                t_final,
                dt
            );

        writeCaseCsv(
            result,
            case_id
        );

        printCaseSummary(
            result,
            case_id
        );
    }

    std::cout
        << "Wrote HW3 validation data to data/validation/hw3/\n";

    return 0;
}
