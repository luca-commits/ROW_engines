#include "eddycurrent.h"
#include "utils.h"
#include <cmath>



int main(int argc, char *argv[]) {
    std::filesystem::path here = __FILE__;
    auto mesh_path = here.remove_filename() / "meshes/cylinder_spherical.msh";
    std::map<int, double> tag_to_current{{1, 0}, {2, 1}, {3, 0}}; // 1 -> air, 2 -> cylinder
    std::map<int, double> tag_to_permeability{{1,1.00000037 * MU_0}, {2,0.999994 * MU_0}, {3, 20 * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 1}, {2, 1}, {3, 1}};
    auto [mesh_p, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_permeability, tag_to_conductivity);
    std::cout << "Hii: " << std::endl;
    std::map<unsigned, std::string> number_to_message {{0, "Success"}, {1, "Numerical Issue"}, {2, "No Convergence"}, {3, "Invalid Input"}};


    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh.NumDofs());
    LF_ASSERT_MSG(N_dofs == mesh_p -> NumEntities(2),
                    " N_dofs must agree with number of nodes");
    std::cout << "Solving BVP with " << N_dofs << " FE dofs" << std::endl;
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
    // Right-hand side vector
    Eigen::VectorXd phi(N_dofs);
    phi.setZero();
    double alpha = 1;
    // std::cout << "phi\n" << phi << std::endl;
    ElemMatProvider elemMatProv (cell_permeability);
    ElemVecProvider elemVecProv (cell_current);
    MassMatProvider massMatProv (cell_conductivity);
    // Cell-oriented assembly
    std::cout << "assebling matrix" << std::endl;
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMatProv, A);
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, massMatProv, M);
    lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);



    // Assembly completed: Convert COO matrix A into CRS format using Eigen's
    // internal conversion routines.

    const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
    const Eigen::SparseMatrix<double> M_crs = M.makeSparse();

    // Create an identity matrix of the same size as A_crs
    Eigen::SparseMatrix<double> I(A_crs.rows(), A_crs.cols());
    I.setIdentity();  // This makes the matrix an identity matrix

    const Eigen::SparseMatrix<double> preconditioner_crs = 0 * A_crs + 1 * M_crs;
    const Eigen::SparseMatrix<double> preconditioner_matrix = 0 * A_crs + 1* M_crs;

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > cholesky; 
    cholesky.compute(preconditioner_crs);

    Eigen::VectorXd sol_vec(N_dofs);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double>> preconditioner;
    preconditioner.compute(preconditioner_crs);

    sol_vec.setZero();

    BiCGstab(preconditioner_matrix, N_dofs, preconditioner, phi, sol_vec, 1e-7, 1000, 5);
    double rel_res  = 0;

    rel_res = (preconditioner_matrix * sol_vec - phi).norm() / phi.norm();
    std::cout << "relative residuum " << rel_res << std::endl;
    auto fe_space_p =
     std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    std::cout << "Num Dofs: " << dofh.NumDofs() << std::endl;
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
            discrete_solution[global_idx];
    };
    vtk_writer.WritePointData("A*", *nodal_data);
    lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space_p, discrete_solution);
    std::cout << "computing B from A*" << std::endl;
    utils::MeshFunctionCurl2DFE mf_curl(mf_grad);
    vtk_writer.WriteCellData("B", mf_curl);
    vtk_writer.WriteCellData("test", *nodal_data);

    return 0;
}