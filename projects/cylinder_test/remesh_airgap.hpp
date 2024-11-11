#ifndef REMESH_AIRGAP_HPP
#define REMESH_AIRGAP_HPP

#include <gmsh.h> 
#include <iostream>
#include <Eigen/Dense>
#include <map> 


void remeshAirgap(std::string airgap_geo_file,  const std::string& output_file, double angle) {
    gmsh::initialize();

    gmsh::open(airgap_geo_file);

    
    gmsh::model::geo::rotate({{0, 139}, {0, 140}, {0, 141}, {0, 142}}, 
                              0,
                              0,
                              0,
                              0,
                              0,
                              1,
                              angle);

    // gmsh::fltk::run();

    gmsh::model::mesh::generate(2);

    // Save and show
    gmsh::write(output_file);
    
    gmsh::finalize();
}


unsigned mergeEverything(const std::string& input_file_stator,
                     const std::string& input_file_rotor,
                     const std::string& input_file_airgap,
                     const std::string& output_file) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    unsigned numStableDof;

    try {

        gmsh::option::setNumber("Geometry.Tolerance", 1e-8);
        gmsh::option::setNumber("Geometry.MatchMeshTolerance", 1e-8);
        
        gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
        gmsh::option::setNumber("Mesh.ScalingFactor", 1.0);

        gmsh::merge(input_file_rotor);
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::removeDuplicateNodes();
        gmsh::model::mesh::removeDuplicateElements();

        gmsh::merge(input_file_stator);
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::removeDuplicateNodes();
        gmsh::model::mesh::removeDuplicateElements();

        std::vector<std::size_t> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

        // Number of nodes is the size of nodeTags
        numStableDof = nodeTags.size();

        gmsh::merge(input_file_airgap);
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::removeDuplicateNodes();
        gmsh::model::mesh::removeDuplicateElements();

        
        // gmsh::fltk::run();
        
        gmsh::write(output_file);
        
    } catch (const std::runtime_error& e) {
        gmsh::logger::write("Error during merge operation: " + std::string(e.what()));
        gmsh::finalize();
        throw;
    }

    return numStableDof;
}

#endif //REMESH_AIRGAP_HPP