// #include "permeability.hpp"

#include "magnetic_material.hpp"
#include <iostream>
#include <fstream>

// int main() {
// //     std::ofstream outfile("magnetic_data.txt");
    
// //     if (!outfile.is_open()) {
// //         std::cerr << "Failed to open output file" << std::endl;
// //         return 1;
// //     }
    
// //     outfile << "B H nu dnu/dB\n";
    
// //     const double mu0 = 4 * M_PI * 1e-7;  // Permeability of free space
    
// //     // Use smaller steps for smoother curve
// //     for (double B = 0; B <= 4e-8; B += 4e-11) {
// //         auto [field, nu, dnu_dB] = H(B);
// //         outfile << B << " "<<  field << " " << nu << " " << dnu_dB<< "\n";
// //     }
    
// //     outfile.close();
//     return 0;
// }

int main(){
    std::ofstream outfile("magnetic_data.txt");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }
    MaterialFactory factory;
    auto iron = factory.Create(3);
    
    outfile << "B H nu dnu/dB\n";
    
    const double mu0 = 4 * M_PI * 1e-7;  // Permeability of free space
    
    // Use smaller steps for smoother curve
    for (double B = 0; B <= 2.16; B += 0.01) {
        double nu = iron -> getReluctivity(B);
        double field = iron -> getReluctivity(B) * B;
        double dnu_dB = iron -> getReluctivityDerivative(B);
        outfile << B << " "<<  field << " " << nu << " " <<dnu_dB<< "\n";
    }


    
    outfile.close();

    return 0; 
}