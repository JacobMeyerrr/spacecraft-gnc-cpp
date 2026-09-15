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
// HOMEWORK 2, PROBLEM 3
//
// XPOP-YPSL commanded spacecraft attitude.
//
// (a) Compute and plot commanded roll/pitch/yaw.
// (b) Compute and plot solar-array gimbal angle.
// (c) Extra credit: compute magnetic field direction in body frame.
//
// IMPORTANT SOURCE CHOICES
//
// The submitted solutions are not all mutually consistent:
//
// - Some use different Earth radii / J2 models.
// - Some use different geomagnetic assumptions.
// - The cleanest mutually consistent attitude formulation uses
//   the same two-body orbit/reference frame as Problem 2.
//
// For this implementation:
//
//   mu = 398600 km^3/s^2
//   R_Earth = 6378 km
//   e = 0.25
//   i = 28.5 deg
//   RAAN = 30 deg
//   arg(perigee) = 45 deg
//
// The Sun model follows the lecture-style approximation used by
// the detailed submitted solution:
//   S(t) = S0 + omega_sun * t
//
// Extra credit uses the 11 deg geomagnetic tilt and 20 deg
// geomagnetic longitude model used by the detailed JT/Kim-style
// course solutions. The original assignment only explicitly states
// |B|=1, so this extra-credit geomagnetic model is an assumption.
// ================================================================


// ---------------------------------------------------------------
// COURSE passive/frame rotation matrices.
//
// These are the same sign convention already used by the HW1
// RotationMatrix class.
//
// They are useful here because Problem 3 constructs a DCM whose
// rows are body-axis components in the reference frame.
// ---------------------------------------------------------------

gnc::math::Matrix3 R1_course(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R << 1.0, 0.0, 0.0,
         0.0, c,   s,
         0.0, -s,  c;

    return R;
}

gnc::math::Matrix3 R2_course(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R << c,   0.0, -s,
         0.0, 1.0, 0.0,
         s,   0.0, c;

    return R;
}

gnc::math::Matrix3 R3_course(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);

    gnc::math::Matrix3 R;

    R << c,  s,   0.0,
         -s, c,   0.0,
         0.0, 0.0, 1.0;

    return R;
}


// ---------------------------------------------------------------
// ACTIVE orbit rotation matrices.
//
// These are used to turn the perifocal state into ECI state.
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
// Kepler equation:
//
//     M = E - e sin(E)
//
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
        double f =
            E
            - e * std::sin(E)
            - M;

        double df =
            1.0
            - e * std::cos(E);

        double correction =
            f / df;

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
    // Problem inputs
    // ============================================================

    constexpr double mu =
        398600.0;

    constexpr double earth_radius =
        6378.0;

    constexpr double perigee_altitude =
        600.0;

    constexpr double e =
        0.25;

    constexpr double inclination =
        28.5 * gnc::constants::d2r;

    constexpr double RAAN =
        30.0 * gnc::constants::d2r;

    constexpr double argument_perigee =
        45.0 * gnc::constants::d2r;


    // ------------------------------------------------------------
    // Orbital parameters
    // ------------------------------------------------------------

    double perigee_radius =
        earth_radius
        + perigee_altitude;

    double a =
        perigee_radius
        / (1.0 - e);

    double p =
        a * (1.0 - e*e);

    double mean_motion =
        std::sqrt(
            mu / (a*a*a)
        );

    double orbital_period =
        2.0 * gnc::constants::pi
        / mean_motion;

    double five_orbit_duration =
        5.0 * orbital_period;


    // ------------------------------------------------------------
    // Perifocal -> ECI
    // ------------------------------------------------------------

    gnc::math::Matrix3 Q_PQW_to_ECI =
        R3_active(RAAN)
        * R1_active(inclination)
        * R3_active(argument_perigee);


    // ------------------------------------------------------------
    // Sun model.
    //
    // July 9, 2026 at 16:00 UTC is approximately 110 days and
    // 16 hours after the March 21 vernal-equinox reference used
    // by the lecture approximation.
    //
    // The detailed submitted solution used:
    //
    //     S0 = 360 * days / 365.25
    // ------------------------------------------------------------

    constexpr double days_since_vernal_equinox =
        110.0 + 16.0 / 24.0;

    double sun_angle_0_deg =
        360.0
        * days_since_vernal_equinox
        / 365.25;

    double sun_angle_rate_deg_s =
        360.0
        / (365.25 * 24.0 * 3600.0);

    constexpr double ecliptic_inclination =
        23.45 * gnc::constants::d2r;


    // ------------------------------------------------------------
    // Five hours between prime meridian crossing and perigee.
    //
    // Earth rotates about 15 deg/hour.
    // ------------------------------------------------------------

    double prime_meridian_angle_0_deg =
        15.0 * 5.0;

    double earth_rotation_rate_deg_s =
        360.0
        / (24.0 * 3600.0);


    // ------------------------------------------------------------
    // Extra-credit geomagnetic model.
    //
    // These 11 deg / 20 deg values are used by the detailed
    // course-style submissions. The assignment itself only
    // specifies |B|=1, so this is explicitly an assumption for
    // the extra-credit model.
    // ------------------------------------------------------------

    constexpr double geomagnetic_tilt =
        11.0 * gnc::constants::d2r;

    constexpr double geomagnetic_longitude =
        20.0 * gnc::constants::d2r;


    // ------------------------------------------------------------
    // Time samples
    // ------------------------------------------------------------

    constexpr int N =
        4001;


    std::ofstream output(
        "data/validation/hw2/problem3_output.csv"
    );

    if(!output.is_open())
    {
        std::cerr
            << "Could not open "
            << "data/validation/hw2/problem3_output.csv\n";

        return 1;
    }

    output << std::setprecision(15);

    output
        << "time_s,orbit_number,"
        << "roll_deg,pitch_deg,yaw_deg,"
        << "delta_deg,"
        << "Bx_body,By_body,Bz_body\n";


    double max_kepler_residual =
        0.0;

    double max_xpop_error =
        0.0;

    double max_ypsl_error =
        0.0;

    double max_solar_alignment_error =
        0.0;

    double min_magnetic_norm =
        1.0e300;

    double max_magnetic_norm_error =
        0.0;


    for(int k = 0;
        k < N;
        ++k)
    {
        // --------------------------------------------------------
        // Time
        // --------------------------------------------------------

        double fraction =
            static_cast<double>(k)
            / static_cast<double>(N - 1);

        double time =
            fraction
            * five_orbit_duration;

        double orbit_number =
            time / orbital_period;


        // --------------------------------------------------------
        // Orbit propagation
        // --------------------------------------------------------

        double M =
            mean_motion * time;

        double E =
            solveKepler(M, e);

        double residual =
            E
            - e * std::sin(E)
            - M;

        max_kepler_residual =
            std::max(
                max_kepler_residual,
                std::abs(residual)
            );


        double denominator =
            1.0
            - e * std::cos(E);

        double cos_nu =
            (
                std::cos(E)
                - e
            )
            / denominator;

        double sin_nu =
            std::sqrt(1.0 - e*e)
            * std::sin(E)
            / denominator;

        double radius =
            a
            * (
                1.0
                - e * std::cos(E)
            );


        gnc::math::Vector3 r_pqw;

        r_pqw <<
            radius * cos_nu,
            radius * sin_nu,
            0.0;


        gnc::math::Vector3 v_pqw;

        double velocity_scale =
            std::sqrt(
                mu / p
            );

        v_pqw <<
            -velocity_scale * sin_nu,
            velocity_scale * (e + cos_nu),
            0.0;


        gnc::math::Vector3 r_eci =
            Q_PQW_to_ECI
            * r_pqw;

        gnc::math::Vector3 v_eci =
            Q_PQW_to_ECI
            * v_pqw;


        // --------------------------------------------------------
        // Reference frame:
        //
        //   x_R = direction of motion
        //   y_R = negative orbit normal
        //   z_R = nadir
        // --------------------------------------------------------

        gnc::math::Vector3 r_hat =
            r_eci / r_eci.norm();

        gnc::math::Vector3 h_hat =
            r_eci.cross(v_eci);

        h_hat =
            h_hat / h_hat.norm();

        gnc::math::Vector3 t_hat =
            h_hat.cross(r_hat);

        t_hat =
            t_hat / t_hat.norm();


        // Matrix whose rows transform ECI components into
        // reference-frame components.
        gnc::math::Matrix3 referenceFromECI;

        referenceFromECI.row(0) =
            t_hat.transpose();

        referenceFromECI.row(1) =
            (-h_hat).transpose();

        referenceFromECI.row(2) =
            (-r_hat).transpose();


        // --------------------------------------------------------
        // Sun direction in ECI.
        // --------------------------------------------------------

        double sun_angle_deg =
            sun_angle_0_deg
            + sun_angle_rate_deg_s * time;

        double sun_angle =
            sun_angle_deg
            * gnc::constants::d2r;

        gnc::math::Vector3 sun_eci;

        sun_eci <<
            std::cos(sun_angle),
            std::sin(sun_angle)
                * std::cos(ecliptic_inclination),
            std::sin(sun_angle)
                * std::sin(ecliptic_inclination);

        sun_eci =
            sun_eci
            / sun_eci.norm();


        // --------------------------------------------------------
        // Express Sun line in reference frame.
        // --------------------------------------------------------

        gnc::math::Vector3 sun_reference =
            referenceFromECI
            * sun_eci;


        // ========================================================
        // Problem 3(a)
        //
        // XPOP:
        //
        //     x_B = y_R
        //
        // For the course 1-2-3 DCM, this gives:
        //
        //     theta = 0
        //     psi   = 90 deg
        //
        // YPSL:
        //
        //     y_B dot sun = 0
        //
        // With theta=0 and psi=90, solve:
        //
        //     phi = atan2(s_R1, s_R3)
        // ========================================================

        double phi =
            std::atan2(
                sun_reference(0),
                sun_reference(2)
            );

        double theta =
            0.0;

        double psi =
            90.0 * gnc::constants::d2r;


        // --------------------------------------------------------
        // Course 1-2-3 commanded DCM:
        //
        //     C_R_B = R1(phi) R2(theta) R3(psi)
        //
        // These course R matrices are the same sign convention
        // used by the HW1 RotationMatrix class.
        // --------------------------------------------------------

        gnc::math::Matrix3 C_reference_to_body =
            R1_course(phi)
            * R2_course(theta)
            * R3_course(psi);


        // --------------------------------------------------------
        // Sun line in body coordinates.
        //
        // Since C_reference_to_body is defined from reference
        // coordinates to body coordinates:
        //
        //     s_B = C_R_B * s_R
        // --------------------------------------------------------

        gnc::math::Vector3 sun_body =
            C_reference_to_body
            * sun_reference;


        // XPOP check:
        //
        // x_B should equal +y_R.
        double xpop_check =
            C_reference_to_body(0,1);

        double xpop_error =
            std::abs(
                xpop_check - 1.0
            );

        max_xpop_error =
            std::max(
                max_xpop_error,
                xpop_error
            );


        // YPSL check:
        //
        // y_B dot Sun should be zero.
        double ypsl_error =
            std::abs(
                sun_body(1)
            );

        max_ypsl_error =
            std::max(
                max_ypsl_error,
                ypsl_error
            );


        // ========================================================
        // Problem 3(b)
        //
        // Solar array rotates around body y-axis.
        //
        // At delta = 0:
        //
        //     n_array = +z_B
        //
        // After rotation:
        //
        //     n_array = [sin(delta), 0, cos(delta)]^T
        //
        // Since YPSL puts the Sun line in the body x-z plane:
        //
        //     delta = atan2(s_Bx, s_Bz)
        // ========================================================

        double delta =
            std::atan2(
                sun_body(0),
                sun_body(2)
            );


        gnc::math::Vector3 array_normal_body;

        array_normal_body <<
            std::sin(delta),
            0.0,
            std::cos(delta);


        double alignment_cosine =
            array_normal_body.dot(
                sun_body
            );

        alignment_cosine =
            clampUnit(
                alignment_cosine
            );

        double pointing_error =
            std::acos(
                alignment_cosine
            );

        max_solar_alignment_error =
            std::max(
                max_solar_alignment_error,
                pointing_error
            );


        // ========================================================
        // Problem 3(c) EXTRA CREDIT
        //
        // Simplified geomagnetic model:
        //
        //     B_M = [0,0,1]^T
        //
        // Geographic -> geomagnetic:
        //
        //     C_G_M = R1(epsilon) R3(Delta)
        //
        // Geographic frame rotates relative to ECI with the
        // Earth rotation rate.
        // ========================================================

        double prime_meridian_angle_deg =
            prime_meridian_angle_0_deg
            + earth_rotation_rate_deg_s * time;

        double prime_meridian_angle =
            prime_meridian_angle_deg
            * gnc::constants::d2r;


        gnc::math::Matrix3 C_G_M =
            R1_course(
                geomagnetic_tilt
            )
            * R3_course(
                geomagnetic_longitude
            );


        gnc::math::Vector3 B_geomagnetic;

        B_geomagnetic <<
            0.0,
            0.0,
            1.0;


        // B_G = C_G_M^T * B_M
        gnc::math::Vector3 B_geographic =
            C_G_M.transpose()
            * B_geomagnetic;


        // The assignment says perigee occurs 5 hours after the
        // prime meridian crossed the inertial vernal-equinox axis.
        //
        // Course frame transformation:
        //
        //     C_G_I = R3(alpha_e)
        //
        // therefore:
        //
        //     B_I = C_G_I^T B_G
        gnc::math::Matrix3 C_G_ECI =
            R3_course(
                prime_meridian_angle
            );

        gnc::math::Vector3 B_eci =
            C_G_ECI.transpose()
            * B_geographic;

        B_eci =
            B_eci
            / B_eci.norm();


        // ECI -> reference -> body.
        gnc::math::Vector3 B_reference =
            referenceFromECI
            * B_eci;

        gnc::math::Vector3 B_body =
            C_reference_to_body
            * B_reference;

        B_body =
            B_body
            / B_body.norm();


        double B_norm =
            B_body.norm();

        min_magnetic_norm =
            std::min(
                min_magnetic_norm,
                B_norm
            );

        max_magnetic_norm_error =
            std::max(
                max_magnetic_norm_error,
                std::abs(B_norm - 1.0)
            );


        // --------------------------------------------------------
        // CSV output.
        // --------------------------------------------------------

        output
            << time << ","
            << orbit_number << ","

            << phi * gnc::constants::r2d << ","
            << theta * gnc::constants::r2d << ","
            << psi * gnc::constants::r2d << ","

            << delta * gnc::constants::r2d << ","

            << B_body(0) << ","
            << B_body(1) << ","
            << B_body(2)
            << "\n";
    }


    output.close();


    // ------------------------------------------------------------
    // Final terminal summary.
    // ------------------------------------------------------------

    std::cout
        << std::setprecision(12);

    std::cout
        << "HW2 Problem 3\n\n";

    std::cout
        << "Semi-major axis: "
        << a
        << " km\n";

    std::cout
        << "Orbital period: "
        << orbital_period
        << " s\n";

    std::cout
        << "Initial Sun angle: "
        << sun_angle_0_deg
        << " deg\n";

    std::cout
        << "Initial prime-meridian angle: "
        << prime_meridian_angle_0_deg
        << " deg\n\n";


    std::cout
        << "Problem 3(a):\n";

    // Recompute the first sample cleanly for display.
    {
        double E0 = 0.0;

        double denom0 =
            1.0 - e;

        double cos_nu0 =
            1.0;

        double sin_nu0 =
            0.0;

        double radius0 =
            a * (1.0 - e * std::cos(E0));

        gnc::math::Vector3 r_pqw0;
        r_pqw0 << radius0, 0.0, 0.0;

        gnc::math::Vector3 v_pqw0;

        double vscale0 =
            std::sqrt(mu / p);

        v_pqw0 <<
            0.0,
            vscale0 * (e + 1.0),
            0.0;

        (void)denom0;
        (void)cos_nu0;
        (void)sin_nu0;

        gnc::math::Vector3 r0 =
            Q_PQW_to_ECI
            * r_pqw0;

        gnc::math::Vector3 v0 =
            Q_PQW_to_ECI
            * v_pqw0;

        gnc::math::Vector3 rhat0 =
            r0 / r0.norm();

        gnc::math::Vector3 hhat0 =
            r0.cross(v0);

        hhat0 =
            hhat0 / hhat0.norm();

        gnc::math::Vector3 that0 =
            hhat0.cross(rhat0);

        that0 =
            that0 / that0.norm();

        gnc::math::Matrix3 referenceFromECI0;

        referenceFromECI0.row(0) =
            that0.transpose();

        referenceFromECI0.row(1) =
            (-hhat0).transpose();

        referenceFromECI0.row(2) =
            (-rhat0).transpose();

        double S0 =
            sun_angle_0_deg
            * gnc::constants::d2r;

        gnc::math::Vector3 sun0;

        sun0 <<
            std::cos(S0),
            std::sin(S0)
                * std::cos(ecliptic_inclination),
            std::sin(S0)
                * std::sin(ecliptic_inclination);

        sun0 =
            sun0 / sun0.norm();

        gnc::math::Vector3 sunR0 =
            referenceFromECI0
            * sun0;

        double phi0 =
            std::atan2(
                sunR0(0),
                sunR0(2)
            );

        std::cout
            << "Roll phi(0) = "
            << phi0 * gnc::constants::r2d
            << " deg\n";

        std::cout
            << "Pitch theta(0) = 0 deg\n";

        std::cout
            << "Yaw psi(0) = 90 deg\n";
    }

    std::cout
        << "\nProblem 3(b):\n";

    // The detailed submitted solution reported approximately
    // 6.07 deg initially and about 6.07-6.18 deg over 5 orbits.
    std::cout
        << "See problem3_output.csv / Python plot for delta(t).\n";

    std::cout
        << "Expected detailed-solution range: approximately "
        << "6.07 to 6.18 deg.\n";

    std::cout
        << "\nProblem 3(c), extra credit:\n";

    std::cout
        << "Minimum |B_B| = "
        << min_magnetic_norm
        << "\n";

    std::cout
        << "Maximum magnetic-field norm error = "
        << max_magnetic_norm_error
        << "\n";


    std::cout
        << "\nValidation checks:\n";

    std::cout
        << "Maximum Kepler residual = "
        << max_kepler_residual
        << " rad\n";

    std::cout
        << "Maximum XPOP error = "
        << max_xpop_error
        << "\n";

    std::cout
        << "Maximum YPSL error = "
        << max_ypsl_error
        << "\n";

    std::cout
        << "Maximum solar-array pointing error = "
        << max_solar_alignment_error
            * gnc::constants::r2d
        << " deg\n";

    std::cout
        << "\nWrote:\n"
        << "data/validation/hw2/problem3_output.csv\n";

    return 0;
}