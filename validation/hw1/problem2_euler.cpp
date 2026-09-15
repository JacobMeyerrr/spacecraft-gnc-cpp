#include "gnc/attitude/RotationMatrix.hpp"
#include "gnc/math/Constants.hpp"

#include <iostream>

int main()
{

    // Homework 1, Problem 2, Part 2a
    double angle1 = 30.0 * gnc::constants::d2r;
    double angle2 = -40.0 * gnc::constants::d2r;
    double angle3 = 75.0 * gnc::constants::d2r;

    gnc::math::Matrix3 R =
        gnc::attitude::RotationMatrix::fromEulerSequence(
            231,
            angle1,
            angle2,
            angle3
        );

    std::cout << "R_231:\n";
    std::cout << R << "\n";

    // Homework 1, Problem 2, Part 2b

    // R_321 = [ c_psi*c_theta,                                  s_psi*c_theta,                                  -s_theta,
    //          -c_phi*s_psi + s_phi*s_theta*c_psi,              c_phi*c_psi + s_phi*s_theta*s_psi,              s_phi*c_theta,
    //           s_phi*s_psi + c_phi*s_theta*c_psi,             -s_phi*c_psi + c_phi*s_theta*s_psi,              c_phi*c_theta ]

   // Therefore:
    // psi   = atan2(R(1,2), R(1,1))
    // theta = asin(-R(1,3))
    // phi   = -atan2(R(2,3), R(3,3))

    // Evaluate
    double psi = std::atan2(R(0,1), R(0,0));
    double theta = std::asin(-R(0,2));
    double phi = -std::atan2(R(1,2), R(2,2));

    // Convert radians to degrees
    psi *= gnc::constants::r2d;
    theta *= gnc::constants::r2d;
    phi *= gnc::constants::r2d;

    // Print Results
    std::cout << "\nEquivalent R_321 Euler Angles:\n";
    std::cout << "Psi (Yaw)   = " << psi << " deg\n";
    std::cout << "Theta (Pitch) = " << theta << " deg\n";
    std::cout << "Phi (Roll)   = " << phi << " deg\n";

    return 0;
}