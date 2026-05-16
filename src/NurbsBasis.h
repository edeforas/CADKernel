#pragma once

#include <vector>
#include <array>

namespace NurbsBasis {

// Mathematical constants
constexpr double Pi = 3.14159265358979323846;

// Numerical tolerances
constexpr double EpsilonSmall = 1.e-20;
constexpr double EpsilonNumerical = 1.e-14;
constexpr double EpsilonGeometric = 1.e-6;

// Default linear knot vector for degree-1 operations
const std::array<double, 4> LinearKnots = {0., 0., 1., 1.};





// Compute basis function derivatives up to 2nd order
// ders should be sized [3][degree+1]
void basis_function_derivatives(
    int degree, 
    const std::vector<double>& knots, 
    int span, 
    double u, 
    std::vector<std::vector<double>>& ders);

// Overload for fixed-size array (for compatibility with existing code)
void basis_function_derivatives(
    int degree, 
    const std::vector<double>& knots, 
    int span, 
    double u, 
    double ders[3][32]);

// Compute only the basis functions (0th derivatives)
void basis_functions(
    int span, 
    double u, 
    int degree, 
    const std::vector<double>& knots, 
    std::vector<double>& N);

} // namespace NurbsBasis