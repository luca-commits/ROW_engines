#ifndef PERMEABILITY_HPP
#define PERMEABILITY_HPP

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

#include <fstream> 

using namespace alglib;

double permeability(double H)
{
    try
    {
        //
        // We use cubic spline to interpolate f(x)=x^2 sampled 
        // at 5 equidistant nodes on [-1,+1].
        //
        // First, we use default boundary conditions ("parabolically terminated
        // spline") because cubic spline built with such boundary conditions 
        // will exactly reproduce any quadratic f(x).
        //
        // Then we try to use natural boundary conditions
        //     d2S(-1)/dx^2 = 0.0
        //     d2S(+1)/dx^2 = 0.0
        // and see that such spline interpolated f(x) with small error.
        //

        std::ifstream infile("tabular_permeability.txt");
        std::string x_str, y_str;

        if (infile.is_open()) {
            std::getline(infile, x_str);
            std::cout << "X string: " << x_str << std::endl;
            std::getline(infile, y_str);
            std::cout << "Y string: " << y_str << std::endl;
            infile.close();
        } else {
            std::cerr << "Error: Unable to open tabular_permeability.txt" << std::endl;
            // You might want to check if the file exists in the correct directory
            std::cerr << "Current working directory might not be what you expect" << std::endl;
        }

        real_1d_array x(x_str.c_str());
        real_1d_array y(y_str.c_str());

        double v;
        spline1dinterpolant s;
        ae_int_t natural_bound_type = 2;
        //
        // Test exact boundary conditions: build S(x), calculare S(0.25)
        // (almost same as original function)
        //
        spline1dbuildcubic(x, y, s);
        v = spline1dcalc(s, H);
        printf("%.4f\n", double(v)); // EXPECTED: 0.0625
        return v; 

    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }
    return 0;
}


#endif //PERMEABILITY_HPP