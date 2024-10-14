#include "eddycurrent.h"
#include "utils.h"
#include <cmath>

#define MU_0 4 * M_PI * 1e-7

int main(int argc, char *argv[]) {
    std::filesystem::path here = __FILE__;
    auto mesh_path = here.remove_filename() / "meshes/cylinder_spherical.msh";
    std::map<int, double> tag_to_current{{1, 0}, {2, 1}, {3, 0}}; // 1 -> air, 2 -> cylinder
    std::map<int, double> tag_to_permeability{{1,1.00000037}, {2,0.999994}, {3, 1}};
    std::map<int, double> tag_to_conductivity{{1, 1}, {2, 1}, {3, 1}};
    auto [mesh_p, cell_current, cell_permeability, cell_conductivity] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_permeability, tag_to_conductivity);
    Eigen::VectorXd discrete_solution = eddycurrent::solve(mesh_p, cell_current, cell_permeability, cell_conductivity);
    auto fe_space_p =
     std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    
    // Concatenate the desired path and filename
    std::string vtk_filename = std::string("vtk_files/eddy_solution_") + argv[1] + std::string(".vtk");

    // Pass the concatenated string to VtkWrit
    lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
    vtk_writer.WriteCellData("current", cell_current);
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < discrete_solution.rows();
        global_idx++) {
        nodal_data->operator()(dofh.Entity(global_idx)) =
            discrete_solution[global_idx];
    };
    vtk_writer.WritePointData("A*", *nodal_data);
    lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space_p, discrete_solution);
    std::cout << "computing B from A*" << std::endl;
    utils::MeshFunctionCurl2DFE mf_curl(mf_grad);
    vtk_writer.WriteCellData("B", mf_curl);

    return 0;
}