#include "gnc/models/EarthPointingSpacecraft.hpp"

#include "gnc/attitude/AttitudeKinematics.hpp"
#include "gnc/attitude/EulerAngles.hpp"
#include "gnc/dynamics/GravityGradient.hpp"
#include "gnc/math/Constants.hpp"

#include <cmath>

namespace gnc::models {

    EarthPointingSpacecraft::EarthPointingSpacecraft(
        double earth_mu,
        double earth_radius,
        double altitude,
        const gnc::math::Matrix3& inertia)
        : mu_(earth_mu),
          radius_(earth_radius + altitude),
          omega0_(std::sqrt(
              earth_mu / std::pow(radius_, 3)
          )),
          orbital_period_(
              gnc::constants::two_pi / omega0_
          ),
          inertia_(inertia),
          D_(omega0_ * (gnc::math::Matrix3() <<
              0.0, -2.0, -10.0,
              2.0, 0.0, -4.0,
              10.0, 4.0, 0.0
          ).finished()),
          K_(omega0_ * omega0_ * (gnc::math::Matrix3() <<
              320.0, -6.0, -1.0,
              -8.0, 150.0, -1.0,
              -4.0, 3.0, 30.0
          ).finished())
    {
        const gnc::math::Vector3 omega_eq(
            0.0,
            -omega0_,
            0.0
        );

        const gnc::math::Vector3 nadir_reference(
            0.0,
            0.0,
            1.0
        );

        const gnc::math::Vector3 gyro_eq =
            omega_eq.cross(inertia_ * omega_eq);

        const gnc::math::Vector3 gg_eq =
            gnc::dynamics::gravityGradientTorque(
                omega0_,
                nadir_reference,
                inertia_
            );

        trim_torque_ =
            gyro_eq - gg_eq;
    }

    double EarthPointingSpacecraft::orbitalRadius() const
    {
        return radius_;
    }

    double EarthPointingSpacecraft::orbitalRate() const
    {
        return omega0_;
    }

    double EarthPointingSpacecraft::orbitalPeriod() const
    {
        return orbital_period_;
    }

    const gnc::math::Matrix3&
    EarthPointingSpacecraft::inertia() const
    {
        return inertia_;
    }

    const gnc::math::Matrix3&
    EarthPointingSpacecraft::linearDampingMatrix() const
    {
        return D_;
    }

    const gnc::math::Matrix3&
    EarthPointingSpacecraft::linearStiffnessMatrix() const
    {
        return K_;
    }

    const gnc::math::Vector3&
    EarthPointingSpacecraft::trimTorque() const
    {
        return trim_torque_;
    }

    gnc::math::Vector6
    EarthPointingSpacecraft::nonlinearDerivative(
        double,
        const gnc::math::Vector6& state,
        const gnc::math::Vector3& disturbance) const
    {
        const gnc::math::Vector3 eta =
            state.head<3>();

        const gnc::math::Vector3 omega =
            state.tail<3>();

        const gnc::math::Matrix3 R =
            gnc::attitude::dcm123({eta(0), eta(1), eta(2)});

        const gnc::math::Vector3 e2(
            0.0,
            1.0,
            0.0
        );

        const gnc::math::Vector3 e3(
            0.0,
            0.0,
            1.0
        );

        const gnc::math::Vector3 omega_body_reference =
            omega + omega0_ * R * e2;

        const gnc::math::Vector3 eta_dot =
            gnc::attitude::differentialKinematics123(eta)
            * omega_body_reference;

        const gnc::math::Vector3 nadir_body =
            R * e3;

        const gnc::math::Vector3 gg =
            gnc::dynamics::gravityGradientTorque(
                omega0_,
                nadir_body,
                inertia_
            );

        const gnc::math::Vector3 omega_dot =
            inertia_.inverse() * (
                trim_torque_
                + disturbance
                + gg
                - omega.cross(inertia_ * omega)
            );

        gnc::math::Vector6 derivative;
        derivative.head<3>() = eta_dot;
        derivative.tail<3>() = omega_dot;
        return derivative;
    }

    gnc::math::Vector6
    EarthPointingSpacecraft::linearDerivative(
        const gnc::math::Vector6& state,
        const gnc::math::Vector3& disturbance) const
    {
        const gnc::math::Vector3 eta =
            state.head<3>();

        const gnc::math::Vector3 eta_dot =
            state.tail<3>();

        const gnc::math::Vector3 eta_ddot =
            inertia_.inverse() * (
                disturbance
                - D_ * eta_dot
                - K_ * eta
            );

        gnc::math::Vector6 derivative;
        derivative.head<3>() = eta_dot;
        derivative.tail<3>() = eta_ddot;
        return derivative;
    }

}
