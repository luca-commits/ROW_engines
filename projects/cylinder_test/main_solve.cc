#include "eddycurrent.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    std::cout << "A ------"<< std::endl;
    std::filesystem::path here = __FILE__;
     auto mesh_path = here.remove_filename() / "meshes/cylinder.msh";

    std::map<int, double> tag_to_current{{1, 0}, {2, 10000000}};
    std::map<int, double> tag_to_permeability{{1,1}, {2,1}};
    auto [mesh_p, cell_current, cell_permeability] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_permeability);
    Eigen::VectorXd discrete_solution = eddycurrent::solve(mesh_p, cell_current, cell_permeability);
    auto fe_space_p =
    std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    
    // Concatenate the desired path and filename
    std::string vtk_filename = std::string("vtk_files/eddy_solution_") + argv[1] + std::string(".vtk");

    // Pass the concatenated string to VtkWriter
    lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
    vtk_writer.WriteCellData("current", cell_current);
    //std::cout << "discrete solution" << discrete_solution << std::endl;
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
    //lf::mesh::utils::CodimMeshDataSet<double>

    vtk_writer.WriteCellData("B", mf_curl);

    // Eigen::Matrix <double, 2, 3> V;
    // V <<  0, 0, 1, 
    //       0, 1, 0; 
    // std::cout << "gradBaryCoords Test\n" << eddycurrent::GradsBaryCoords(V) << std::endl;
    // Eigen::Matrix <double , 2 , 3> X = eddycurrent::GradsBaryCoords(V);
    // Eigen::Matrix <double, 2, 3> temp = X;
    // temp.block<1, 3>(0, 0) = X.block<1,3> (1, 0);
    // temp.block<1, 3>(1, 0) = X.block<1,3> (0, 0);
    // double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
    // // std::cout << "element matrix\n" << X.transpose() * temp << std::endl;
    // std::cout << area * temp.transpose() * temp; 
    // return 0;
}