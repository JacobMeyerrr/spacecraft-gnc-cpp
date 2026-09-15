#include "gnc/simulation/Simulator.hpp"
#include "gnc/simulation/Integrator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gnc::simulation {

    SimulationHistory simulateFixedStep(
        double t0,
        double tf,
        double dt,
        const gnc::math::Vector6& x0,
        const std::function<gnc::math::Vector6(
            double,
            const gnc::math::Vector6&)>& derivative)
    {
        if(dt <= 0.0 || tf < t0)
        {
            throw std::invalid_argument(
                "Invalid simulation time interval"
            );
        }

        SimulationHistory result;
        gnc::math::Vector6 state = x0;
        double time = t0;

        result.time.push_back(time);
        result.state.push_back(state);

        while(time < tf)
        {
            const double h =
                std::min(dt, tf - time);

            state = rk4Step(
                derivative,
                time,
                state,
                h
            );

            time += h;

            result.time.push_back(time);
            result.state.push_back(state);
        }

        return result;
    }

}
