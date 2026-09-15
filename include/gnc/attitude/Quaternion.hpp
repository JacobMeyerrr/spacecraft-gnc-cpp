#pragma once

#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::attitude {

    class Quaternion
    {
        public:
            // Identity quaternion.
            //
            // q = [0 0 0 1]^T
            //
            // The course convention is vector-first / scalar-last.
            Quaternion();

            // Direct constructor.
            //
            // q1, q2, q3 = vector/imaginary part
            // q4        = scalar/real part
            Quaternion(
                double q1,
                double q2,
                double q3,
                double q4
            );

            // Create a quaternion from an axis-angle rotation.
            //
            // q_vec = n_hat * sin(theta/2)
            // q4    = cos(theta/2)
            static Quaternion fromAxisAngle(
                const gnc::math::Vector3& axis,
                double angle
            );

            // Convert a 3x3 rotation matrix into a quaternion.
            static Quaternion fromRotationMatrix(
                const gnc::math::Matrix3& R
            );

            // Convert the COURSE Euler-sequence rotation matrix
            // into a quaternion.
            //
            // This is intentionally done through the course
            // RotationMatrix class because the basic R1/R2/R3
            // matrices in this project use the course frame/
            // direction-cosine convention.
            static Quaternion fromEulerSequence(
                int sequence,
                double angle1,
                double angle2,
                double angle3
            );

            // Convert quaternion to the standard active rotation
            // matrix used for rotating vectors:
            //
            //     v_new = R * v_old
            gnc::math::Matrix3 toRotationMatrix() const;

            // Rotate a vector with:
            //
            //     v_rotated = q * [v,0] * q*
            //
            // for a unit quaternion.
            gnc::math::Vector3 rotate(
                const gnc::math::Vector3& v
            ) const;

            // Quaternion normalization.
            void normalize();

            Quaternion normalized() const;

            // Quaternion conjugate:
            //
            // q* = [-q_vec, q4]
            Quaternion conjugate() const;

            // General quaternion inverse:
            //
            // q^-1 = q* / ||q||^2
            Quaternion inverse() const;

            // Axis-angle representation.
            double angle() const;
            gnc::math::Vector3 axis() const;

            // Hamilton quaternion product.
            Quaternion operator*(
                const Quaternion& other
            ) const;

            // Accessors.
            double q1() const;
            double q2() const;
            double q3() const;
            double q4() const;

        private:
            double q_1;
            double q_2;
            double q_3;
            double q_4;
    };

}