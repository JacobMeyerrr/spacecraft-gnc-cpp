#include "gnc/attitude/Quaternion.hpp"
#include "gnc/attitude/RotationMatrix.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

int main()
{
    std::cout << std::setprecision(8);

    // ============================================================
    // HOMEWORK 2, PROBLEM 1
    // ============================================================

    // ------------------------------------------------------------
    // Part 1(a)
    //
    // Compute a quaternion for a 25 degree rotation about:
    //
    //     a = 2 i - 3 j + 1 k
    // ------------------------------------------------------------

    gnc::math::Vector3 axis_a;

    axis_a <<
        2.0,
        -3.0,
        1.0;

    gnc::attitude::Quaternion q_a =
        gnc::attitude::Quaternion::fromAxisAngle(
            axis_a,
            25.0 * gnc::constants::d2r
        );

    std::cout
        << "P1(a) Quaternion [q1 q2 q3 q4]^T:\n"
        << q_a.q1() << "\n"
        << q_a.q2() << "\n"
        << q_a.q3() << "\n"
        << q_a.q4() << "\n\n";

    // Some submitted solutions print scalar-first.
    // That is ONLY a display convention.
    std::cout
        << "P1(a) Same quaternion in scalar-first display [q4 q1 q2 q3]^T:\n"
        << q_a.q4() << "\n"
        << q_a.q1() << "\n"
        << q_a.q2() << "\n"
        << q_a.q3() << "\n\n";


    // ------------------------------------------------------------
    // Part 1(b)
    //
    // Given:
    //
    // q = [-0.3388, 0.11290, 0.2259, 0.9063]^T
    // ------------------------------------------------------------

    gnc::attitude::Quaternion q_b(
        -0.3388,
        0.11290,
        0.2259,
        0.9063
    );

    double angle_b =
        q_b.angle();

    gnc::math::Vector3 axis_b =
        q_b.axis();

    std::cout
        << "P1(b) Angle [deg]: "
        << angle_b * gnc::constants::r2d
        << "\n";

    std::cout
        << "P1(b) Axis:\n"
        << axis_b
        << "\n\n";


    // ------------------------------------------------------------
    // Part 1(c)
    //
    // Course Euler sequence:
    //
    //     R_312(35 deg, -50 deg, 70 deg)
    //
    //     R_312 = R2(70) R1(-50) R3(35)
    // ------------------------------------------------------------

    gnc::math::Matrix3 R_c =
        gnc::attitude::RotationMatrix::fromEulerSequence(
            312,
            35.0 * gnc::constants::d2r,
            -50.0 * gnc::constants::d2r,
            70.0 * gnc::constants::d2r
        );

    gnc::attitude::Quaternion q_c =
        gnc::attitude::Quaternion::fromRotationMatrix(
            R_c
        );

    double angle_c =
        q_c.angle();

    gnc::math::Vector3 axis_c =
        q_c.axis();

    std::cout
        << "P1(c) R_312 rotation matrix:\n"
        << R_c
        << "\n\n";

    std::cout
        << "P1(c) Equivalent angle [deg]: "
        << angle_c * gnc::constants::r2d
        << "\n";

    std::cout
        << "P1(c) Equivalent axis:\n"
        << axis_c
        << "\n\n";


    // ------------------------------------------------------------
    // Part 1(d)
    //
    // b = [2,-3,1]^T
    //
    // Three successive ACTIVE axis-angle rotations:
    //
    // 1) +35 deg about [-2,4,1]
    // 2) -65 deg about [2,5,-2]
    // 3) +120 deg about [-5,4,3]
    // ------------------------------------------------------------

    gnc::math::Vector3 b;

    b <<
        2.0,
        -3.0,
        1.0;

    gnc::math::Vector3 axis1;

    axis1 <<
        -2.0,
        4.0,
        1.0;

    gnc::math::Vector3 axis2;

    axis2 <<
        2.0,
        5.0,
        -2.0;

    gnc::math::Vector3 axis3;

    axis3 <<
        -5.0,
        4.0,
        3.0;


    gnc::attitude::Quaternion q1 =
        gnc::attitude::Quaternion::fromAxisAngle(
            axis1,
            35.0 * gnc::constants::d2r
        );

    gnc::attitude::Quaternion q2 =
        gnc::attitude::Quaternion::fromAxisAngle(
            axis2,
            -65.0 * gnc::constants::d2r
        );

    gnc::attitude::Quaternion q3 =
        gnc::attitude::Quaternion::fromAxisAngle(
            axis3,
            120.0 * gnc::constants::d2r
        );


    // Rotation 1 happens first.
    // Rotation 3 happens last.
    //
    // Therefore:
    //
    //     q_total = q3 * q2 * q1

    gnc::attitude::Quaternion q_total =
        q3 * q2 * q1;


    // ------------------------------------------------------------
    // Part 1(d)(A)
    // ------------------------------------------------------------

    gnc::math::Vector3 b_new =
        q_total.rotate(b);

    std::cout
        << "P1(d)(A) Rotated vector:\n"
        << b_new
        << "\n\n";


    // ------------------------------------------------------------
    // Part 1(d)(B)
    // ------------------------------------------------------------

    double angle_d =
        q_total.angle();

    gnc::math::Vector3 axis_d =
        q_total.axis();

    std::cout
        << "P1(d)(B) Equivalent angle [deg]: "
        << angle_d * gnc::constants::r2d
        << "\n";

    std::cout
        << "P1(d)(B) Equivalent axis:\n"
        << axis_d
        << "\n\n";


    // ------------------------------------------------------------
    // Part 1(d)(C)
    // ------------------------------------------------------------

    gnc::math::Matrix3 R_d =
        q_total.toRotationMatrix();

    std::cout
        << "P1(d)(C) Resulting rotation matrix:\n"
        << R_d
        << "\n\n";


    // ------------------------------------------------------------
    // Part 1(d)(D)
    //
    // Equivalent 321 Euler angles for the ACTIVE matrix:
    //
    //     yaw   = atan2(R21,R11)
    //     pitch = asin(-R13)
    //     roll  = -atan2(R23,R33)
    // ------------------------------------------------------------

    double yaw =
        std::atan2(
            R_d(0,1),
            R_d(0,0)
        );

    double pitch =
        std::asin(
            std::max(
                -1.0,
                std::min(
                    1.0,
                    -R_d(0,2)
                )
            )
        );

    double roll =
        -std::atan2(
            R_d(1,2),
            R_d(2,2)
        );

    std::cout
        << "P1(d)(D) Equivalent 321 Euler angles [deg]:\n";

    std::cout
        << "Roll (phi)   = "
        << roll * gnc::constants::r2d
        << "\n";

    std::cout
        << "Pitch (theta) = "
        << pitch * gnc::constants::r2d
        << "\n";

    std::cout
        << "Yaw (psi)     = "
        << yaw * gnc::constants::r2d
        << "\n";

    return 0;
}