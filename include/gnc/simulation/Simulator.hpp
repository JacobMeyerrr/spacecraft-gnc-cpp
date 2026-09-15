#pragma once

#include "gnc/math/Vector.hpp"

#include <functional>
#include <vector>

namespace gnc::simulation {

    struct SimulationHistory
    {
        std::vector<double> time;
        std::vector<gnc::math::Vector6> state;
    };

    SimulationHistory simulateFixedStep(
        double t0,
        double tf,
        double dt,
        const gnc::math::Vector6& x0,
        const std::function<gnc::math::Vector6(
            double,
            const gnc::math::Vector6&)>& derivative);

}
