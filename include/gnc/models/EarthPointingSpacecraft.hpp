#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::models {

    class EarthPointingSpacecraft
    {
        public:
            EarthPointingSpacecraft(
                double earth_mu,
                double earth_radius,
                double altitude,
                const gnc::math::Matrix3& inertia);

            double orbitalRadius() const;
            double orbitalRate() const;
            double orbitalPeriod() const;

            const gnc::math::Matrix3& inertia() const;
            const gnc::math::Matrix3& linearDampingMatrix() const;
            const gnc::math::Matrix3& linearStiffnessMatrix() const;
            const gnc::math::Vector3& trimTorque() const;

            gnc::math::Vector6 nonlinearDerivative(
                double t,
                const gnc::math::Vector6& state,
                const gnc::math::Vector3& disturbance) const;

            gnc::math::Vector6 linearDerivative(
                const gnc::math::Vector6& state,
                const gnc::math::Vector3& disturbance) const;

        private:
            double mu_;
            double radius_;
            double omega0_;
            double orbital_period_;
            gnc::math::Matrix3 inertia_;
            gnc::math::Matrix3 D_;
            gnc::math::Matrix3 K_;
            gnc::math::Vector3 trim_torque_;
    };

}
