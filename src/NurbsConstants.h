#pragma once

#include <array>
#include <vector>

namespace NurbsConstants {

// Mathematical constants
constexpr double Pi = 3.14159265358979323846;

// Numerical tolerances
constexpr double EpsilonSmall = 1.e-20;
constexpr double EpsilonNumerical = 1.e-14;
constexpr double EpsilonGeometric = 1.e-6;

// Default linear knot vector for degree-1 operations
const std::array<double, 4> LinearKnots = {0., 0., 1., 1.};

} // namespace NurbsConstants