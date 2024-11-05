#ifndef REMESH_AIRGAP_HPP
#define REMESH_AIRGAP_HPP

#include <gmsh.h> 
#include <iostream>
#include <Eigen/Dense>
#include <map> 

void remeshAirgap(const std::string& input_file_stator, const std::string& input_file_rotor, const std::string& output_file, double angle) {
    gmsh::initialize();

    // Load rotor mesh and get nodes
    gmsh::open("airgap.geo");
    
    gmsh::model::geo::rotate({{0, 39}, {0, 40}, {0, 41}, {0, 42}}, 
                              0,
                              0,
                              0,
                              0,
                              0,
                              1,
                              angle);

    gmsh::model::mesh::generate(2);

    gmsh::fltk::run();
    // Save and show
    gmsh::write(output_file);
    
    gmsh::finalize();
}


void mergeEverything(const std::string& input_file_stator, 
                    const std::string& input_file_rotor, 
                    const std::string& input_file_airgap, 
                    const std::string& output_file) {
    gmsh::initialize();
    
    try {
        // Set binary format for better precision
        gmsh::option::setNumber("Mesh.Binary", 1);
        
        // Store original entities before merge
        std::vector<std::pair<int, int>> dimTags;
        
        // Merge files - gmsh::merge() throws an exception if it fails
        gmsh::open(input_file_rotor);
        // Get and store entity tags after first merge
        gmsh::model::getEntities(dimTags);
        gmsh::option::setNumber("Geometry.Tolerance", 1e-7);
        gmsh::merge(input_file_stator);
        gmsh::merge(input_file_airgap);
        
        // Ensure model coherence
        gmsh::model::mesh::removeDuplicateNodes();
        
        // Optional: Unify mesh if needed
        gmsh::option::setNumber("Mesh.Binary", 1);
        gmsh::option::setNumber("Mesh.SaveAll", 1);
        
        // Display the merged mesh
        gmsh::fltk::run();
        
        // Write with specific format
        gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
        gmsh::write(output_file);
        
    } catch (const std::runtime_error& e) {
        gmsh::logger::write("Error during merge operation: " + std::string(e.what()));
        gmsh::finalize();
        throw;
    }
    
    gmsh::finalize();
}





#endif //REMESH_AIRGAP_HPP