#include "NumericalOptimization.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

double testfunction(double x1, double x2)
{
    return 5 * x1 + pow(x1, 2) + 3 * x2 + pow(x2, 2);
}
int main() {
    auto local_function = [](std::vector<double> input)
    {
        return testfunction(input[0], input[1]);
    };
    std::vector<double> init = {0, 1};
    std::vector<double> out { numerical::nelder_mead_fmin(local_function, init) };
    std::cout << "function = 5 * x1 + pow(x1, 2) + 3 * x2 + pow(x2, 2)" << std::endl;
    std::cout << "x1 " << out[0] << std::endl;
    std::cout << "x2 " << out[1] << std::endl;
    std::cout << "min " << testfunction(out[0], out[1]) << std::endl;
}
