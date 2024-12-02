#ifndef MAGNETIC_MATERIAL_HPP
#define MAGNETIC_MATERIAL_HPP

#include <cmath>
#include <memory>
#include <map>
#include <set>
#include <stdexcept>
#include <lf/assemble/assemble.h>

#include "interpolation.h"

// Abstract base class for magnetic material behavior using reluctivity
class MagneticMaterial {
public:
    virtual ~MagneticMaterial() = default;
    // Returns H given B (inverse of original B(H) function)
    // Returns reluctivity (ν = H/B) at given B
    virtual Eigen::Vector2d getH(const Eigen::Vector2d& B) const = 0; 
    virtual double getReluctivity(const double B) const = 0;
    virtual double getReluctivityDerivative(double B) const = 0;
protected:
    static constexpr double mu0_ = 4 * M_PI * 1e-7; // Vacuum permeability
    static constexpr double nu0_ = 1.0 / mu0_;      // Vacuum reluctivity
};

// Linear material (constant reluctivity)
class LinearMaterial : public MagneticMaterial {
public:
    // Constructor now takes relative reluctivity (nu_r = 1/mu_r)
    LinearMaterial(double nu_r = 1.0) : nu_r_(nu_r) {}

    Eigen::Vector2d getH(const Eigen::Vector2d& B) const override {
        return nu0_ * nu_r_ * B;
    }

    double getReluctivity(const double /*B*/) const override {
        return nu0_ * nu_r_;
    }

     double getReluctivityDerivative(double B_magnitude) const override{
        return 0; 
    }

private:
    const double nu_r_;  
};

class FerromagneticMaterial : public MagneticMaterial {
public:
    FerromagneticMaterial(){
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
                    std::getline(infile, y_str);
                    infile.close();
                } else {
                    std::cerr << "Error: Unable to open tabular_permeability.txt" << std::endl;
                    // You might want to check if the file exists in the correct directory
                    std::cerr << "Current working directory might not be what you expect" << std::endl;
                }


                alglib::real_1d_array x(x_str.c_str());
                alglib::real_1d_array y(y_str.c_str());

                // Build cubic spline with specified boundary condition
                alglib::spline1dbuildlinear(x, y, s_);

                
            }
            catch(alglib::ap_error alglib_exception)
            {
                printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
                throw std::runtime_error("error in construction of spline");
            }
    }

    Eigen::Vector2d getH(const Eigen::Vector2d& B) const override {
        return getReluctivity(B.norm()) * B; 
    }


    double getReluctivity(double B) const override {
        if (B < 1e-12) {
            double reluctivity = alglib::spline1dcalc(s_, 1e-10) / 1e-10; 
            // std::cout << "reluctivity: " << reluctivity << std::endl; 
            return reluctivity;  // Return vacuum reluctivity for very small B
        }

        double dx = 0;
        double function = 0;
        double dx2 = 0;
        alglib::spline1ddiff(s_, B, function, dx, dx2);
        // std::cout << "B : " << B << std::endl; 
        // std::cout << "reluctivity : " << alglib::spline1dcalc(s_, B) / B << std::endl;
        return function / B;
    }
    
    double getReluctivityDerivative(double B) const override{

        if (B < 1e-12) {
            return 0.0;
            // throw std::runtime_error("B is smaller than ")
        }


        double dx = 0;
        double function = 0;
        double dx2 = 0;
        alglib::spline1ddiff(s_, B, function, dx, dx2);
        
        // std::cout << "reluctivity derivative : " << 1e10  * ((B / (B * B + 1)) - std::atan(B)) / (B * B) << std::endl;  
        // return ((beta_ * B / (beta_ * beta_ * B * B + 1)) - std::atan(beta_ * B)) / (B * B);
        // std::cout << "B : " << B << std::endl; 
        // std::cout << "dx";
        // std::cout << "reluctivity derivative : " <<  dx/B - function/(B*B) << std::endl; 
        // if (dx/B - function/(B*B) < 0) throw std::runtime_error("derivative smaller than zero");
        return dx/B - function/(B*B); 
    }

private:
    alglib::spline1dinterpolant s_;
};

// Factory to create materials
class MaterialFactory {
public:
    static std::shared_ptr<MagneticMaterial> Create(int material_tag) {

        switch(material_tag) {
            case 3: // ring
                // throw std::runtime_error("found a ring element");
                return std::make_shared<FerromagneticMaterial>();
            case 2: // cylinder
                return std::make_shared<LinearMaterial>(1/(0.999994));  // nu_r = 1/mu_r
            case 1: //air
                return std::make_shared<LinearMaterial>(1/(1.00000037)); 
            case 4: //airgap
                return std::make_shared<LinearMaterial>(1/(1.00000037));                 
            default:
                throw std::runtime_error("Unknown material tag: " + std::to_string(material_tag));
        }
    }
};

#endif //MAGNETIC_MATERIAL_HPP