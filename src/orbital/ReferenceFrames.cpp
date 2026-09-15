#include "gnc/orbital/ReferenceFrames.hpp"

namespace gnc::orbital {

    ReferenceFrame earthPointingFrame(
        const gnc::math::Vector3& position_eci,
        const gnc::math::Vector3& velocity_eci)
    {
        ReferenceFrame frame;

        const gnc::math::Vector3 r_hat =
            position_eci.normalized();

        const gnc::math::Vector3 h_hat =
            position_eci.cross(velocity_eci).normalized();

        const gnc::math::Vector3 t_hat =
            h_hat.cross(r_hat).normalized();

        frame.x_eci = t_hat;
        frame.y_eci = -h_hat;
        frame.z_eci = -r_hat;

        frame.eci_to_reference.row(0) =
            frame.x_eci.transpose();

        frame.eci_to_reference.row(1) =
            frame.y_eci.transpose();

        frame.eci_to_reference.row(2) =
            frame.z_eci.transpose();

        return frame;
    }

}
