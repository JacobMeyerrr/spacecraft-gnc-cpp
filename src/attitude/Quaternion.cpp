#include "gnc/attitude/Quaternion.hpp"
#include "gnc/attitude/RotationMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gnc::attitude {

    Quaternion::Quaternion()
        : q_1(0.0),
          q_2(0.0),
          q_3(0.0),
          q_4(1.0)
    {
    }


    Quaternion::Quaternion(
        double q1,
        double q2,
        double q3,
        double q4)
        : q_1(q1),
          q_2(q2),
          q_3(q3),
          q_4(q4)
    {
    }


    Quaternion Quaternion::fromAxisAngle(
        const gnc::math::Vector3& axis,
        double angle)
    {
        // Always normalize the supplied axis first.
        double axis_norm = axis.norm();

        if(axis_norm == 0.0)
        {
            throw std::invalid_argument(
                "Rotation axis cannot have zero magnitude"
            );
        }

        gnc::math::Vector3 n_hat =
            axis / axis_norm;

        double half_angle =
            angle / 2.0;

        double s =
            std::sin(half_angle);

        double c =
            std::cos(half_angle);

        // Lecture convention:
        //
        // q = [ n_hat*sin(theta/2)
        //       cos(theta/2) ]
        //
        // i.e. q1-q3 are the vector part and q4 is scalar.
        return Quaternion(
            n_hat(0) * s,
            n_hat(1) * s,
            n_hat(2) * s,
            c
        );
    }


    Quaternion Quaternion::fromRotationMatrix(
        const gnc::math::Matrix3& R)
    {
        // Standard active-vector quaternion/DCM relationship.

        double q4_squared =
            0.25 * (
                1.0
                + R.trace()
            );

        if(q4_squared > 1.0e-12)
        {
            double q4 =
                std::sqrt(q4_squared);

            double q1 =
                (R(2,1) - R(1,2))
                / (4.0 * q4);

            double q2 =
                (R(0,2) - R(2,0))
                / (4.0 * q4);

            double q3 =
                (R(1,0) - R(0,1))
                / (4.0 * q4);

            return Quaternion(
                q1,
                q2,
                q3,
                q4
            ).normalized();
        }

        // Near 180 degrees, use whichever vector component is
        // largest to avoid dividing by a tiny number.

        double q1_squared =
            0.25 * (
                1.0
                + R(0,0)
                - R(1,1)
                - R(2,2)
            );

        double q2_squared =
            0.25 * (
                1.0
                - R(0,0)
                + R(1,1)
                - R(2,2)
            );

        double q3_squared =
            0.25 * (
                1.0
                - R(0,0)
                - R(1,1)
                + R(2,2)
            );

        if(q1_squared >= q2_squared &&
           q1_squared >= q3_squared &&
           q1_squared > 1.0e-12)
        {
            double q1 =
                std::sqrt(q1_squared);

            double q2 =
                (R(0,1) + R(1,0))
                / (4.0 * q1);

            double q3 =
                (R(0,2) + R(2,0))
                / (4.0 * q1);

            double q4 =
                (R(2,1) - R(1,2))
                / (4.0 * q1);

            return Quaternion(
                q1,
                q2,
                q3,
                q4
            ).normalized();
        }

        if(q2_squared >= q1_squared &&
           q2_squared >= q3_squared &&
           q2_squared > 1.0e-12)
        {
            double q2 =
                std::sqrt(q2_squared);

            double q1 =
                (R(0,1) + R(1,0))
                / (4.0 * q2);

            double q3 =
                (R(1,2) + R(2,1))
                / (4.0 * q2);

            double q4 =
                (R(0,2) - R(2,0))
                / (4.0 * q2);

            return Quaternion(
                q1,
                q2,
                q3,
                q4
            ).normalized();
        }

        if(q3_squared > 1.0e-12)
        {
            double q3 =
                std::sqrt(q3_squared);

            double q1 =
                (R(0,2) + R(2,0))
                / (4.0 * q3);

            double q2 =
                (R(1,2) + R(2,1))
                / (4.0 * q3);

            double q4 =
                (R(1,0) - R(0,1))
                / (4.0 * q3);

            return Quaternion(
                q1,
                q2,
                q3,
                q4
            ).normalized();
        }

        throw std::invalid_argument(
            "Rotation matrix cannot be converted to quaternion"
        );
    }


    Quaternion Quaternion::fromEulerSequence(
        int sequence,
        double angle1,
        double angle2,
        double angle3)
    {
        // IMPORTANT:
        //
        // Do not directly multiply elementary quaternions here
        // using the active-vector convention.
        //
        // The course Euler-sequence matrices use the course frame /
        // DCM convention. For example:
        //
        //     R_312 = R2(angle3) * R1(angle2) * R3(angle1)
        //
        // Therefore build that exact course matrix first and then
        // convert the matrix to a quaternion.

        gnc::math::Matrix3 R =
            RotationMatrix::fromEulerSequence(
                sequence,
                angle1,
                angle2,
                angle3
            );

        return fromRotationMatrix(R);
    }


    gnc::math::Matrix3 Quaternion::toRotationMatrix() const
    {
        // Normalize first so numerical roundoff does not violate
        // the unit-quaternion condition.

        Quaternion q =
            normalized();

        double q1 = q.q_1;
        double q2 = q.q_2;
        double q3 = q.q_3;
        double q4 = q.q_4;

        gnc::math::Matrix3 R;

        // Standard ACTIVE rotation matrix:
        //
        //     v_new = R * v_old

        R << q4*q4 + q1*q1 - q2*q2 - q3*q3,
             2.0*(q1*q2 - q3*q4),
             2.0*(q1*q3 + q2*q4),

             2.0*(q1*q2 + q3*q4),
             q4*q4 - q1*q1 + q2*q2 - q3*q3,
             2.0*(q2*q3 - q1*q4),

             2.0*(q1*q3 - q2*q4),
             2.0*(q2*q3 + q1*q4),
             q4*q4 - q1*q1 - q2*q2 + q3*q3;

        return R;
    }


    gnc::math::Vector3 Quaternion::rotate(
        const gnc::math::Vector3& v) const
    {
        Quaternion q =
            normalized();

        Quaternion v_q(
            v(0),
            v(1),
            v(2),
            0.0
        );

        // Quaternion rotation:
        //
        //     v' = q * v * q*

        Quaternion result =
            q
            * v_q
            * q.conjugate();

        gnc::math::Vector3 rotated;

        rotated <<
            result.q_1,
            result.q_2,
            result.q_3;

        return rotated;
    }


    void Quaternion::normalize()
    {
        double norm =
            std::sqrt(
                q_1*q_1
                + q_2*q_2
                + q_3*q_3
                + q_4*q_4
            );

        if(norm == 0.0)
        {
            throw std::invalid_argument(
                "Cannot normalize zero quaternion"
            );
        }

        q_1 /= norm;
        q_2 /= norm;
        q_3 /= norm;
        q_4 /= norm;
    }


    Quaternion Quaternion::normalized() const
    {
        Quaternion q =
            *this;

        q.normalize();

        return q;
    }


    Quaternion Quaternion::conjugate() const
    {
        return Quaternion(
            -q_1,
            -q_2,
            -q_3,
            q_4
        );
    }


    Quaternion Quaternion::inverse() const
    {
        double norm_squared =
            q_1*q_1
            + q_2*q_2
            + q_3*q_3
            + q_4*q_4;

        if(norm_squared == 0.0)
        {
            throw std::invalid_argument(
                "Cannot invert zero quaternion"
            );
        }

        Quaternion q =
            conjugate();

        q.q_1 /= norm_squared;
        q.q_2 /= norm_squared;
        q.q_3 /= norm_squared;
        q.q_4 /= norm_squared;

        return q;
    }


    double Quaternion::angle() const
    {
        Quaternion q =
            normalized();

        double q4 =
            std::clamp(
                q.q_4,
                -1.0,
                1.0
            );

        return 2.0 * std::acos(q4);
    }


    gnc::math::Vector3 Quaternion::axis() const
    {
        Quaternion q =
            normalized();

        double sin_half_angle =
            std::sqrt(
                std::max(
                    0.0,
                    1.0 - q.q_4*q.q_4
                )
            );

        gnc::math::Vector3 axis;

        if(sin_half_angle < 1.0e-12)
        {
            axis <<
                0.0,
                0.0,
                0.0;

            return axis;
        }

        axis <<
            q.q_1 / sin_half_angle,
            q.q_2 / sin_half_angle,
            q.q_3 / sin_half_angle;

        return axis;
    }


    Quaternion Quaternion::operator*(
        const Quaternion& other) const
    {
        // Hamilton product.
        //
        // q = [x y z w]^T
        // where w = q4.

        double x1 = q_1;
        double y1 = q_2;
        double z1 = q_3;
        double w1 = q_4;

        double x2 = other.q_1;
        double y2 = other.q_2;
        double z2 = other.q_3;
        double w2 = other.q_4;

        return Quaternion(
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2,
            w1*w2 - x1*x2 - y1*y2 - z1*z2
        );
    }


    double Quaternion::q1() const
    {
        return q_1;
    }


    double Quaternion::q2() const
    {
        return q_2;
    }


    double Quaternion::q3() const
    {
        return q_3;
    }


    double Quaternion::q4() const
    {
        return q_4;
    }

}