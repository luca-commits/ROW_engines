#ifndef PERMEABILITY_HPP
#define PERMEABILITY_HPP

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

#include <fstream> 

using namespace alglib;

// double H(double B)
// {
//     try
//     {
//         //
//         // We use cubic spline to interpolate f(x)=x^2 sampled 
//         // at 5 equidistant nodes on [-1,+1].
//         //
//         // First, we use default boundary conditions ("parabolically terminated
//         // spline") because cubic spline built with such boundary conditions 
//         // will exactly reproduce any quadratic f(x).
//         //
//         // Then we try to use natural boundary conditions
//         //     d2S(-1)/dx^2 = 0.0
//         //     d2S(+1)/dx^2 = 0.0
//         // and see that such spline interpolated f(x) with small error.
//         //
//         std::cout << "B " << B << std::endl; 
//         std::ifstream infile("tabular_permeability.txt");
//         std::string x_str, y_str;

//         if (infile.is_open()) {
//             std::getline(infile, x_str);
//             std::cout << "x string: " << std::endl << x_str << std::endl;
//             std::getline(infile, y_str);
//             std::cout << "y string: " << std::endl  << y_str << std::endl;
//             infile.close();
//         } else {
//             std::cerr << "Error: Unable to open tabular_permeability.txt" << std::endl;
//             // You might want to check if the file exists in the correct directory
//             std::cerr << "Current working directory might not be what you expect" << std::endl;
//         }


//         real_1d_array x(x_str.c_str());
//         real_1d_array y(y_str.c_str());
//         double v;
//         spline1dinterpolant s;

//         // Set boundary conditions
//         // 2 means we're setting derivative value
//         // 0.0 is the derivative value we want (zero in this case)
//         ae_int_t bound_type_left = 1;  // specify derivative at left
//         ae_int_t bound_type_right = 0; // natural spline at right
//         double bound_l = 0.0;          // derivative = 0 at left boundary
//         double bound_r = 0.0;          // not used for natural boundary condition

//         // Build cubic spline with specified boundary condition
//         spline1dbuildmonotone(x, y, s);

//         v = spline1dcalc(s, B);

//         std::cout << "V " << v << std::endl;

//         // if (v < 0) throw std::runtime_error("reluctivity should be positive"); 
//         return v; 

//     }
//     catch(alglib::ap_error alglib_exception)
//     {
//         printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
//         return 1;
//     }
//     return 0;
// }

std::tuple<double, double, double> H(double B)
{
    std::ifstream infile("tabular_permeability.txt");
    std::string x_str, y_str;

    if (infile.is_open()) {
        std::getline(infile, x_str);
        std::getline(infile, y_str);
        infile.close();
    } else {
        std::cerr << "Error: Unable to open tabular_permeability.txt" << std::endl;
        std::cerr << "Current working directory might not be what you expect" << std::endl;
    }

    real_1d_array x(x_str.c_str());
    real_1d_array y(y_str.c_str());
    spline1dinterpolant s;

    // Build linear interpolant
    spline1dbuildlinear(x, y, s);

    double dx = 0;
    double function = 0;
    double dx2 = 0;
    alglib::spline1ddiff(s, B, function, dx, dx2);

    return {function, function / B, dx/B - function / (B * B)};
}
#endif //PERMEABILITY_HPP