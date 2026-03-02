#ifndef LinearAlgebra_
#define LinearAlgebra_

#include <vector>

bool solve_linear_system(std::vector<double> M, const std::vector<double>& b, std::vector<double>& x, int n);

#endif