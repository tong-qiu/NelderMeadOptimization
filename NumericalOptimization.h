/*
numerical optimization lib

Version 1.0
Release data: 27 Mar 2018
Copyright 2018 Tong Qiu 

The program contains the following numerical optimization functions:
1. nelder_mead() 
   The function use Nelder-Mead simplex algorithm to find the min or max
   value of a function

*/

#ifndef NUMERICAL_OPTIMIZATION_H
#define NUMERICAL_OPTIMIZATION_H

#include<string>
#include<vector>
#include<functional>
namespace numerical
{
    // the function uses Nelder-Mead simplex algorithm to find the minimum or maximum
    // value of a function
    // return the parameter at the min or max value
    // input_function: a point to a function whose parameter is a vector
    // start: inital guess of best point
    // method: find min or max
    // scale: the length of each side of the initial simplex
    // precision: the precision of the value of the min of the function
    // max_step: max iterations
    std::vector<double> nelder_mead_fmin(std::function <double(std::vector<double>)>,
        std::vector <double>, std::string = "min",
        double = 1., double = 1e-4, int = 100000);
}
#endif