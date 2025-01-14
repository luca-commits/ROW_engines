#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"

#include <fstream>

int main (int argc, char *argv[]){

    std::filesystem::path here = __FILE__;
    std::string type_of_rotor = "box_irregular_ring";
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor + "/";

    std::string final_mesh = mesh_path + "rotated_rotor_" + std::to_string(3) + ".msh";
    std::cout << "final mesh : " << final_mesh << std::endl;  
    std::string vtk_filename = std::string("vtk_files/time_dependent/non-linear/eddy_solution_transient_non_linear") + "_" + std::to_string(3) + std::string(".vtk");

    gmsh::initialize();
    gmsh::open(final_mesh); 
    gmsh::fltk::run();
    gmsh::finalize();

    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), final_mesh);
    
    std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};

    lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);

    lf::mesh::utils::CodimMeshDataSet<double> relative_permeability{mesh_p, 0, -1};
    for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
        relative_permeability(*cell) = 1;
    }

    vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
    return 0; 
}