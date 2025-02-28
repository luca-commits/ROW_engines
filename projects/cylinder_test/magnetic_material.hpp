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
        // if (nu_r_ == 0)std::cout << "ring reluctivity : " << nu_r_ << std::endl; 
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
    FerromagneticMaterial() : 
    max_(5000),
    c_(1),
    mu0_(4 * M_PI * 1e-7) {}

    Eigen::Vector2d getH(const Eigen::Vector2d& B) const override {
        return getReluctivity(B.norm()) * B;     
    }


    double getReluctivity(double B) const override {
        if (B < 1e-12) {
            getReluctivity(1e-12);
        }
        double relative_permeability = max_ / (1 + std::pow(B, 4) * max_ / c_) + 1;
        return (1 / relative_permeability) / mu0_;
    }
    
    double getReluctivityDerivative(double B) const override{
        if (B < 1e-12) {
            return 0.0;
        }
        else{
            return (4 * max_ * max_ * c_ * std::pow(B, 3)) / std::pow((max_ * c_ + std::pow(B, 4) * max_ + c_), 2) / mu0_;
        }
    }

private:
    double max_ = 5000; 
    double c_ = 100; 
    double mu0_ = 4 * M_PI * 1e-7;
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
                return std::make_shared<LinearMaterial>(1/(1.));  // nu_r = 1/mu_r
            case 1: //air
                return std::make_shared<LinearMaterial>(1/(1.)); 
            case 4: //airgap
                return std::make_shared<LinearMaterial>(1/(1.));                              
            default:
                throw std::runtime_error("Unknown material tag: " + std::to_string(material_tag));
        }
    }
};

#endif //MAGNETIC_MATERIAL_HPP