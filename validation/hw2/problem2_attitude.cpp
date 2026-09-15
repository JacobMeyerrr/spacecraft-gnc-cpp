#include "gnc/attitude/Quaternion.hpp"
#include "gnc/math/Constants.hpp"
#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// ================================================================
// HOMEWORK 2, PROBLEM 2
//
// Earth-pointing spacecraft attitude relative to ECI.
//
// The submitted solutions that agree with each other use:
//
//   x_R = direction of motion
//   y_R = negative orbit normal
//   z_R = nadir
//
// The fixed body alignment errors are:
//
//   roll  = 20 deg
//   pitch = 15 deg
//   yaw   = 30 deg
//
// We compute the full attitude over 5 orbital periods and write
// the time history to a CSV so Python can make the plots.
// ================================================================


// ---------------------------------------------------------------
// ACTIVE rotation matrices.
//
// These are used for the physical orbit transformation and for
// the fixed attitude error matrix in the body/reference chain.
// ---------------------------------------------------------------

gnc::math::Matrix3 R1_active(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R << 1.0, 0.0, 0.0,
         0.0, c,  -s,
         0.0, s,   c;

    return R;
}

gnc::math::Matrix3 R2_active(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R <<  c, 0.0, s,
         0.0, 1.0, 0.0,
         -s, 0.0, c;

    return R;
}

gnc::math::Matrix3 R3_active(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R << c,  -s, 0.0,
         s,   c, 0.0,
         0.0, 0.0, 1.0;

    return R;
}


// ---------------------------------------------------------------
// Kepler solver:
//
//     M = E - e sin(E)
//
// Newton-Raphson.
// ---------------------------------------------------------------

double solveKepler(
    double M,
    double e)
{
    double E = M;

    for(int iteration = 0;
        iteration < 20;
        ++iteration)
    {
        double residual =
            E
            - e * std::sin(E)
            - M;

        double derivative =
            1.0
            - e * std::cos(E);

        double correction =
            residual / derivative;

        E -= correction;

        if(std::abs(correction) < 1.0e-13)
        {
            break;
        }
    }

    return E;
}


double clampUnit(double value)
{
    return std::max(
        -1.0,
        std::min(
            1.0,
            value
        )
    );
}


int main()
{
    // ============================================================
    // Constants and problem inputs
    // ============================================================

    constexpr double mu =
        398600.0; // km^3/s^2

    constexpr double earth_radius =
        6378.0; // km

    constexpr double perigee_altitude =
        600.0; // km

    constexpr double eccentricity =
        0.25;

    constexpr double inclination =
        28.5 * gnc::constants::d2r;

    constexpr double RAAN =
        30.0 * gnc::constants::d2r;

    constexpr double argument_perigee =
        45.0 * gnc::constants::d2r;

    constexpr double roll_error =
        20.0 * gnc::constants::d2r;

    constexpr double pitch_error =
        15.0 * gnc::constants::d2r;

    constexpr double yaw_error =
        30.0 * gnc::constants::d2r;


    // ------------------------------------------------------------
    // Orbital parameters
    // ------------------------------------------------------------

    double perigee_radius =
        earth_radius
        + perigee_altitude;

    double semi_major_axis =
        perigee_radius
        / (1.0 - eccentricity);

    double semi_latus_rectum =
        semi_major_axis
        * (1.0 - eccentricity * eccentricity);

    double apogee_radius =
        semi_major_axis
        * (1.0 + eccentricity);

    double apogee_altitude =
        apogee_radius
        - earth_radius;

    double mean_motion =
        std::sqrt(
            mu
            / (
                semi_major_axis
                * semi_major_axis
                * semi_major_axis
            )
        );

    double orbital_period =
        2.0 * gnc::constants::pi
        / mean_motion;

    double five_orbit_duration =
        5.0 * orbital_period;


    // ------------------------------------------------------------
    // Perifocal -> ECI active rotation.
    //
    // This follows the standard active vector transformation:
    //
    //     Q = R3(RAAN) R1(i) R3(omega_p)
    // ------------------------------------------------------------

    gnc::math::Matrix3 Q_PQW_to_ECI =
        R3_active(RAAN)
        * R1_active(inclination)
        * R3_active(argument_perigee);


    // ------------------------------------------------------------
    // Fixed alignment-error matrix.
    //
    // This is applied after the nominal Earth-pointing reference
    // frame is built.
    // ------------------------------------------------------------

    gnc::math::Matrix3 C_reference_to_body =
        R3_active(yaw_error)
        * R2_active(pitch_error)
        * R1_active(roll_error);


    // ------------------------------------------------------------
    // Use the same sampling density as the detailed submitted
    // MATLAB solution: 4001 points over exactly 5 periods.
    // ------------------------------------------------------------

    constexpr int number_of_points =
        4001;


    std::ofstream output(
        "data/validation/hw2/problem2_output.csv"
    );

    if(!output.is_open())
    {
        std::cerr
            << "Could not open "
            << "data/validation/hw2/problem2_output.csv\n";

        return 1;
    }

    output
        << std::setprecision(15);

    output
        << "time_s,orbit_number,"
        << "x_eci_km,y_eci_km,z_eci_km,"
        << "vx_eci_km_s,vy_eci_km_s,vz_eci_km_s,"
        << "roll_deg,pitch_deg,yaw_deg,"
        << "q1,q2,q3,q4,"
        << "body_x_eci_1,body_x_eci_2,body_x_eci_3,"
        << "body_y_eci_1,body_y_eci_2,body_y_eci_3,"
        << "body_z_eci_1,body_z_eci_2,body_z_eci_3\n";


    double max_kepler_residual = 0.0;
    double max_orthogonality_error = 0.0;
    double min_determinant = 1.0e300;
    double max_quaternion_matrix_error = 0.0;

    gnc::math::Matrix3 initial_attitude =
        gnc::math::Matrix3::Identity();

    gnc::math::Matrix3 final_attitude =
        gnc::math::Matrix3::Identity();

    bool have_previous_quaternion =
        false;

    gnc::attitude::Quaternion previous_q;


    for(int k = 0;
        k < number_of_points;
        ++k)
    {
        // --------------------------------------------------------
        // Time
        // --------------------------------------------------------

        double fraction =
            static_cast<double>(k)
            / static_cast<double>(number_of_points - 1);

        double time =
            fraction
            * five_orbit_duration;

        double orbit_number =
            time / orbital_period;


        // --------------------------------------------------------
        // Mean anomaly and eccentric anomaly
        // --------------------------------------------------------

        double M =
            mean_motion * time;

        double E =
            solveKepler(
                M,
                eccentricity
            );

        double residual =
            E
            - eccentricity * std::sin(E)
            - M;

        max_kepler_residual =
            std::max(
                max_kepler_residual,
                std::abs(residual)
            );


        // --------------------------------------------------------
        // True anomaly.
        //
        // We compute sine and cosine separately so there is no
        // artificial wrap jump at +/- pi.
        // --------------------------------------------------------

        double denominator =
            1.0
            - eccentricity * std::cos(E);

        double cos_nu =
            (
                std::cos(E)
                - eccentricity
            )
            / denominator;

        double sin_nu =
            std::sqrt(
                1.0
                - eccentricity * eccentricity
            )
            * std::sin(E)
            / denominator;


        // --------------------------------------------------------
        // Radius
        // --------------------------------------------------------

        double radius =
            semi_major_axis
            * (
                1.0
                - eccentricity * std::cos(E)
            );


        // --------------------------------------------------------
        // PQW position/velocity
        // --------------------------------------------------------

        gnc::math::Vector3 r_pqw;

        r_pqw <<
            radius * cos_nu,
            radius * sin_nu,
            0.0;


        gnc::math::Vector3 v_pqw;

        double velocity_scale =
            std::sqrt(
                mu / semi_latus_rectum
            );

        v_pqw <<
            -velocity_scale * sin_nu,
            velocity_scale
                * (
                    eccentricity
                    + cos_nu
                ),
            0.0;


        // --------------------------------------------------------
        // ECI state
        // --------------------------------------------------------

        gnc::math::Vector3 r_eci =
            Q_PQW_to_ECI
            * r_pqw;

        gnc::math::Vector3 v_eci =
            Q_PQW_to_ECI
            * v_pqw;


        // --------------------------------------------------------
        // Build the nominal Earth-pointing reference frame:
        //
        //   x_R = direction of motion
        //   y_R = negative orbit normal
        //   z_R = nadir
        //
        // These three columns are the reference axes expressed
        // in ECI coordinates.
        // --------------------------------------------------------

        gnc::math::Vector3 r_hat =
            r_eci / r_eci.norm();

        gnc::math::Vector3 h_vec =
            r_eci.cross(v_eci);

        gnc::math::Vector3 h_hat =
            h_vec / h_vec.norm();

        gnc::math::Vector3 t_hat =
            h_hat.cross(r_hat);

        t_hat =
            t_hat / t_hat.norm();


        gnc::math::Matrix3 C_ECI_reference;

        C_ECI_reference.col(0) =
            t_hat;

        C_ECI_reference.col(1) =
            -h_hat;

        C_ECI_reference.col(2) =
            -r_hat;


        // --------------------------------------------------------
        // Actual body attitude:
        //
        //     C_ECI_B =
        //         C_ECI_reference * C_reference_to_body
        //
        // Columns are x_B, y_B, z_B expressed in ECI.
        // --------------------------------------------------------

        gnc::math::Matrix3 C_ECI_body =
            C_ECI_reference
            * C_reference_to_body;


        // --------------------------------------------------------
        // Euler-angle extraction.
        //
        // For:
        //
        // C = R3(yaw) R2(pitch) R1(roll)
        //
        // using the active-vector matrices above:
        //
        // pitch = asin(-C31)
        // yaw   = atan2(C21,C11)
        // roll  = atan2(C32,C33)
        // --------------------------------------------------------

        double pitch =
            std::asin(
                clampUnit(
                    -C_ECI_body(2,0)
                )
            );

        double yaw =
            std::atan2(
                C_ECI_body(1,0),
                C_ECI_body(0,0)
            );

        double roll =
            std::atan2(
                C_ECI_body(2,1),
                C_ECI_body(2,2)
            );


        // --------------------------------------------------------
        // Quaternion representation of the body attitude.
        // --------------------------------------------------------

        gnc::attitude::Quaternion q =
            gnc::attitude::Quaternion::fromRotationMatrix(
                C_ECI_body
            );

        // q and -q represent the same attitude.
        // Flip the sign between neighboring samples if necessary
        // so the plotted history is continuous.
        if(have_previous_quaternion)
        {
            double dot_product =
                q.q1()*previous_q.q1()
                + q.q2()*previous_q.q2()
                + q.q3()*previous_q.q3()
                + q.q4()*previous_q.q4();

            if(dot_product < 0.0)
            {
                q =
                    gnc::attitude::Quaternion(
                        -q.q1(),
                        -q.q2(),
                        -q.q3(),
                        -q.q4()
                    );
            }
        }

        previous_q =
            q;

        have_previous_quaternion =
            true;


        // --------------------------------------------------------
        // Quaternion -> matrix round-trip check.
        // --------------------------------------------------------

        gnc::math::Matrix3
            C_from_quaternion =
                q.toRotationMatrix();

        double quaternion_matrix_error =
            (
                C_ECI_body
                - C_from_quaternion
            ).norm();

        max_quaternion_matrix_error =
            std::max(
                max_quaternion_matrix_error,
                quaternion_matrix_error
            );


        // --------------------------------------------------------
        // Rotation matrix health checks.
        // --------------------------------------------------------

        double orthogonality_error =
            (
                C_ECI_body.transpose()
                * C_ECI_body
                - gnc::math::Matrix3::Identity()
            ).norm();

        double determinant =
            C_ECI_body.determinant();

        max_orthogonality_error =
            std::max(
                max_orthogonality_error,
                orthogonality_error
            );

        min_determinant =
            std::min(
                min_determinant,
                determinant
            );


        if(k == 0)
        {
            initial_attitude =
                C_ECI_body;
        }

        if(k == number_of_points - 1)
        {
            final_attitude =
                C_ECI_body;
        }


        // --------------------------------------------------------
        // Save one row to CSV.
        // --------------------------------------------------------

        output
            << time << ","
            << orbit_number << ","

            << r_eci(0) << ","
            << r_eci(1) << ","
            << r_eci(2) << ","

            << v_eci(0) << ","
            << v_eci(1) << ","
            << v_eci(2) << ","

            << roll * gnc::constants::r2d << ","
            << pitch * gnc::constants::r2d << ","
            << yaw * gnc::constants::r2d << ","

            << q.q1() << ","
            << q.q2() << ","
            << q.q3() << ","
            << q.q4() << ",";

        // Body x axis
        output
            << C_ECI_body(0,0) << ","
            << C_ECI_body(1,0) << ","
            << C_ECI_body(2,0) << ",";

        // Body y axis
        output
            << C_ECI_body(0,1) << ","
            << C_ECI_body(1,1) << ","
            << C_ECI_body(2,1) << ",";

        // Body z axis
        output
            << C_ECI_body(0,2) << ","
            << C_ECI_body(1,2) << ","
            << C_ECI_body(2,2)
            << "\n";
    }


    output.close();


    // ------------------------------------------------------------
    // Exact five-orbit closure check.
    // ------------------------------------------------------------

    double closure_error =
        (
            final_attitude
            - initial_attitude
        ).norm();


    std::cout
        << std::setprecision(12);

    std::cout
        << "HW2 Problem 2\n\n";

    std::cout
        << "Semi-major axis: "
        << semi_major_axis
        << " km\n";

    std::cout
        << "Orbital period: "
        << orbital_period
        << " s\n";

    std::cout
        << "Apogee altitude: "
        << apogee_altitude
        << " km\n\n";

    std::cout
        << "Initial C_ECI_B:\n"
        << initial_attitude
        << "\n\n";


    double initial_pitch =
        std::asin(
            clampUnit(
                -initial_attitude(2,0)
            )
        );

    double initial_yaw =
        std::atan2(
            initial_attitude(1,0),
            initial_attitude(0,0)
        );

    double initial_roll =
        std::atan2(
            initial_attitude(2,1),
            initial_attitude(2,2)
        );

    std::cout
        << "Initial Euler angles:\n";

    std::cout
        << "Roll  = "
        << initial_roll * gnc::constants::r2d
        << " deg\n";

    std::cout
        << "Pitch = "
        << initial_pitch * gnc::constants::r2d
        << " deg\n";

    std::cout
        << "Yaw   = "
        << initial_yaw * gnc::constants::r2d
        << " deg\n\n";


    std::cout
        << "Maximum Kepler residual = "
        << max_kepler_residual
        << " rad\n";

    std::cout
        << "Maximum orthogonality error = "
        << max_orthogonality_error
        << "\n";

    std::cout
        << "Minimum determinant = "
        << min_determinant
        << "\n";

    std::cout
        << "Maximum quaternion/matrix round-trip error = "
        << max_quaternion_matrix_error
        << "\n";

    std::cout
        << "Five-orbit closure error = "
        << closure_error
        << "\n\n";

    std::cout
        << "Wrote:\n"
        << "data/validation/hw2/problem2_output.csv\n";

    return 0;
}