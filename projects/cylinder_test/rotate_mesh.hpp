#ifndef ROTATE_MESH_HPP
#define ROTATE_MESH_HPP

#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <gmsh.h>
#include <iostream> 

void rotateAllNodes_alt(const std::string& input_file, const std::string& output_file, double angle) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    
    // Calculate rotation matrix coefficients once
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);
    
    // Open the mesh file
    gmsh::merge(input_file);

    // Get all nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;
    
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);
    
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        // Get original coordinates
        double x = coords[i * 3];
        double y = coords[i * 3 + 1];
        double z = coords[i * 3 + 2];
        
        // Calculate rotated coordinates
        double new_x = x * cosA - y * sinA;
        double new_y = x * sinA + y * cosA;
        
        // Update node coordinates using setNode
        gmsh::model::mesh::setNode(
            nodeTags[i],
            {new_x, new_y, z} ,             
            parametricCoords             
        );
    }

    gmsh::model::mesh::rebuildNodeCache(false);
    gmsh::write(output_file);

    

// Finalize GMSH
gmsh::finalize();
}

void deform_airgap(const std::string& input_file, const std::string& output_file, double angle, double r, double R) {
    // angle = 0.01;
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    
    // Calculate rotation matrix coefficients once
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);
    
    // Open the mesh file
    gmsh::merge(input_file);

    // Get all nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;
    
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);


    
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        // std::cout << "i " << i << std::endl;
        // Get original coordinates
        double x = coords[i * 3];
        double y = coords[i * 3 + 1];
        double z = coords[i * 3 + 2];

        //rotation angle depends on distance from outer border
        double r_x = std::sqrt(x * x + y * y);
        double phi_x = ((r_x - r) / (R - r)) * angle; 
        // std::cout << "r_x : " << r_x << std::endl; 
        // std::cout << "phi_x : " << phi_x << std::endl; 
        double cosA = std::cos(phi_x);
        double sinA = std::sin(phi_x);
        
        // Calculate rotated coordinates
        double new_x = x * cosA - y * sinA;
        double new_y = x * sinA + y * cosA;
        
        // Update node coordinates using setNode
        gmsh::model::mesh::setNode(
            nodeTags[i],
            {new_x, new_y, z} ,             
            parametricCoords             
        );
    }


    gmsh::model::mesh::rebuildNodeCache(false);
    // gmsh::fltk::run();
    gmsh::write(output_file);

    
    // Finalize GMSH
    gmsh::finalize();
}


#endif // ROTATE_MESH_HPP
