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

#include "NumericalOptimization.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
using namespace std;

// the function use Nelder-Mead simplex algorithm to find the min or max
// value of a function
// return the parameter at the min or max value
// input_function: a point to a function whose parameter is a vector
// start: inital guess of best point
// method: find min or max
// scale: the scale of the function
// precision: the precision of the value of the min of the function
// max_step: max iterations

//vector<double> numerical::nelder_mead_fmin(double(*input_function)(vector<double>),
vector<double> numerical::nelder_mead_fmin(function <double(vector<double>)> input_function,
    vector <double> start, string method,
    double scale, double precision, int max_step)
{
    // check the method parameter
    if (method != "min" && method != "max")
    {
        cout << "Error: the method of nelder_mead function has to be min or max" << endl;
        exit(1);
    }

    // the program is designed to find the minimum value.
    // minus sign is added if user wishes to find max value, 
    auto function = [=](double &result, vector<double> parameter)
    {
        if (method == "min")
            result = input_function(parameter);
        else
            result = -input_function(parameter);
    };

    // the total number of the parameters
    int n;
    n = start.size();
    // all vertices of the simplex
    vector<vector<double>> vertices;
    // values of the funtion at each simplex
    double *values = new double[n + 1];

    // first vertex of the simplex
    vertices.push_back(start);
    // calculate the value of the function at the first vertex
    function(values[0], start);
    // construct the rest of the simplex
    for (int i{ 0 }; i < n; i++)
    {
        vector<double> tem_vertices{ start };
        tem_vertices[i] = tem_vertices[i] + scale;
        vertices.push_back(tem_vertices);
        // calculate the value of the function at the vertex
        function(values[i + 1], tem_vertices);
    }

    // find the max, second max and min values and their locations
    double max_value{ *max_element(values, values + n + 1) };
    double min_value{ *min_element(values, values + n + 1) };
    int max_location{ int(distance(values, max_element(values, values + n + 1)) )};
    int min_location{ int(distance(values, min_element(values, values + n + 1)) )};
    double second_max_value{ 0 };
    int second_max_location{ 0 };
    for (int i{ 0 }; i < n + 1; i++)
    {
        if (max_location == i)
            continue;
        if (values[i] >= second_max_value)
        {
            second_max_value = values[i];
            second_max_location = i;
        }
    }

    // variable to count steps
    int step{ 1 };
    // begin optimization
    while (max_value - min_value > precision && step != max_step)
    {
        step++;
        // compute the value at each vertex
        delete[] values;
        values = new double[n + 1];
        for (int i{ 0 }; i < n + 1; i++)
            function(values[i], vertices[i]);
        // find the max, second max and min values and their locations
        max_value = *max_element(values, values + n + 1);
        min_value = *min_element(values, values + n + 1);
        max_location = distance(values, max_element(values, values + n + 1));
        min_location = distance(values, min_element(values, values + n + 1));
        second_max_value = second_max_location = 0;
        for (int i{ 0 }; i < n + 1; i++)
        {
            if (i == max_location)
                continue;
            if (values[i] >= second_max_value)
            {
                second_max_value = values[i];
                second_max_location = i;
            }
        }

        // the worst vertex
        vector<double> worst_vertex{ vertices[max_location] };
        // the best vertex
        vector<double> best_vertex{ vertices[min_location] };

        // calculate the average vertex except the worst
        vector<double> average_vertex;
        for (int i{ 0 }; i < n; i++)
            average_vertex.push_back(0.);
        // loop over vertices
        for (int i{ 0 }; i < n + 1; i++)
        {
            if (i == max_location)
                continue;
            // loop over each coordinate
            for (int j{ 0 }; j < n; j++)
                average_vertex[j] += vertices[i][j] / double(n);
        }

        // reflection
        // calculate reflected vertex and its value
        vector<double> reflected_vertex;
        double reflected_value;
        for (int i{ 0 }; i < n; i++)
            reflected_vertex.push_back(average_vertex[i] +
            (average_vertex[i] - worst_vertex[i]));
        function(reflected_value, reflected_vertex);

        // if the reflected vertex is better than the best vertex
        // expansion
        if (reflected_value < min_value)
        {
            // calculate expanded vertex and its value
            vector<double> expanded_vertex;
            double expanded_value;
            for (int i{ 0 }; i < n; i++)
                expanded_vertex.push_back(reflected_vertex[i] +
                (reflected_vertex[i] - average_vertex[i]));
            function(expanded_value, expanded_vertex);

            // compare the expansion and reflection
            // replace the worst vertex by the best
            if (expanded_value < reflected_value)
                vertices[max_location] = expanded_vertex;
            else
                vertices[max_location] = reflected_vertex;
            continue;
        }

        // if the reflected vertex is worse than the worst vertex
        if (reflected_value > max_value)
        {
            // inside contraction
            vector<double> inside_contraction_vertex;
            double inside_contraction_value;
            for (int i{ 0 }; i < n; i++)
                inside_contraction_vertex.push_back(average_vertex[i] -
                    0.5 * (average_vertex[i] - worst_vertex[i]));
            function(inside_contraction_value, inside_contraction_vertex);

            // if the inside contraction vertex is better than the worst vertex
            if (inside_contraction_value < max_value)
                // replace the worst vetice by inside contracion vertex
                vertices[max_location] = inside_contraction_vertex;
            // if not
            else
            {
                // shrink
                // loop over vertices
                for (int i{ 0 }; i < n + 1; i++)
                {
                    if (i == min_location)
                        continue;
                    // loop over each coordinate
                    for (int j{ 0 }; j < n; j++)
                        vertices[i][j] = vertices[min_location][j] +
                        0.5 * (vertices[i][j] - vertices[min_location][j]);
                }
            }
            continue;
        }

        // if the reflected vertex is better than the second worst vertex
        if (reflected_value < second_max_value)
        {
            // replace the worst by the reflectd vertex
            vertices[max_location] = reflected_vertex;
            continue;
        }

        // if the reflected vertex is worse than the second worst vertex
        // outside contraction
        vector<double> out_contracted_vertex;
        double out_contracted_value;
        for (int i{ 0 }; i < n; i++)
            out_contracted_vertex.push_back(average_vertex[i] +
                0.5 * (average_vertex[i] - worst_vertex[i]));
        function(out_contracted_value, out_contracted_vertex);
        // if the outside constracted vertece is better than reflected vertex
        if (out_contracted_value < reflected_value)
            vertices[max_location] = out_contracted_vertex;
        // if not
        else
        {
            // shrink
            // loop over vertices
            for (int i{ 0 }; i < n + 1; i++)
            {
                if (i == min_location)
                    continue;
                // loop over each coordinate
                for (int j{ 0 }; j < n; j++)
                    vertices[i][j] = vertices[min_location][j] +
                    0.5 * (vertices[i][j] - vertices[min_location][j]);
            }
        }
    }

    if (step == max_step)
        cout << "Warning: max step reached at numerical::nelder_mead_fmin!"<<endl;

    delete[] values;
    return vertices[min_location];
}
