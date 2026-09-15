#pragma once

namespace gnc::simulation {

    template <typename State, typename RHS>
    State rk4Step(
        const RHS& rhs,
        double t,
        const State& x,
        double dt)
    {
        const State k1 = rhs(t, x);
        const State k2 = rhs(t + 0.5 * dt, x + 0.5 * dt * k1);
        const State k3 = rhs(t + 0.5 * dt, x + 0.5 * dt * k2);
        const State k4 = rhs(t + dt, x + dt * k3);

        return x + (dt / 6.0) *
            (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

}
