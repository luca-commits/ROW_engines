#ifndef MAGNETIC_MATERIAL_HPP
#define MAGNETIC_MATERIAL_HPP

#include <cmath>
#include <memory>
#include <map>
#include <set>
#include <stdexcept>
#include <lf/assemble/assemble.h>

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
    FerromagneticMaterial(double mu_initial) : mu_initial_(mu_initial * mu0_) {}

    Eigen::Vector2d getH(const Eigen::Vector2d& B) const override {
        double B_magnitude = B.norm();
        if (B_magnitude < 1e-10) {
            return B / mu_initial_;
        }
        return 0 * B; //getReluctivity(B_magnitude) * B;
    }


    double getReluctivity(double B) const override {
        if (B < 1e-10) {
            return 1/mu_initial_;
        }
        return std::atan(beta_ * B) / B ;
    }
    
    double getReluctivityDerivative(double B) const override{
        if (B < 1e-10) {
                   return 0.0;
            // throw std::runtime_error("B is smaller than ")
        }
        // std::cout << "reluctivity derivative : " << 1e10  * ((B / (B * B + 1)) - std::atan(B)) / (B * B) << std::endl;  
        return ((beta_ * B / (beta_ * beta_ * B * B + 1)) - std::atan(beta_ * B)) / (B * B);
    }

private:
    const double mu_initial_;
    double beta_ = 1e9; 

};

// Factory to create materials
class MaterialFactory {
public:
    static std::shared_ptr<MagneticMaterial> Create(int material_tag) {

        switch(material_tag) {
            case 3: // ring
                // throw std::runtime_error("found a ring element");
                return std::make_shared<FerromagneticMaterial>(1/(0.999994));
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