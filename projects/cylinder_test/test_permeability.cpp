#include "permeability.hpp"
#include <iostream>
#include <fstream>

int main() {
    std::ofstream outfile("magnetic_data.txt");
    
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }
    
    outfile << "H B mu_r\n";
    
    const double mu0 = 4 * M_PI * 1e-7;  // Permeability of free space
    
    // Use smaller steps for smoother curve
    for (double H = 100; H <= 2000; H += 5) {
        double mu_r = permeability(H);
        double B = mu_r * mu0 * H;  // Calculate B from μr
        
        outfile << H << " " << B << " " << mu_r << "\n";
    }
    
    outfile.close();
    return 0;
}