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
     // Initialize GMSH
    gmsh::initialize();
    
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
            nodeTags[i],                    // node tag
            {new_x, new_y, z} , // new coordinates
            parametricCoords             
        );
    }

    gmsh::model::mesh::rebuildNodeCache(false);
    gmsh::write(output_file);
    

// Finalize GMSH
gmsh::finalize();
}

#endif // ROTATE_MESH_HPP
