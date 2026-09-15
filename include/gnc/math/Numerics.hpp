#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace gnc::math {

    inline double clamp(
        double value,
        double lower,
        double upper)
    {
        return std::max(lower, std::min(upper, value));
    }

    inline double clampUnit(double value)
    {
        return clamp(value, -1.0, 1.0);
    }

    inline bool nearlyEqual(
        double a,
        double b,
        double tolerance = 1.0e-12)
    {
        return std::abs(a - b) <= tolerance;
    }

    inline bool isFinite(double value)
    {
        return std::isfinite(value);
    }

}
