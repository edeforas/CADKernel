#pragma once

#include <vector>

namespace NurbsKnots {

// Normalize knots to [0,1] range
void normalize_to_01(std::vector<double>& knots);

// Scale knots to [0,1] range (same as normalize)
inline void scale_knots(std::vector<double>& knots) {
    normalize_to_01(knots);
}

// Build uniform knots
std::vector<double> build_uniform_knots(int degree, int nbPoints);

// Build clamped uniform knots
std::vector<double> build_clamped_uniform_knots(int degree, int nbPoints);

// Build open uniform knots
std::vector<double> build_open_uniform_knots(int degree, int nbPoints);

// Build segmented quadratic knots
std::vector<double> build_segmented_quadratic_knots(int nbSegments);

} // namespace NurbsKnots