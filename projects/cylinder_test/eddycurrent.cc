#include "eddycurrent.h"
#include "BICG_stab.hpp"
#include <iomanip> // For std::setw
#include <cstdlib> // For rand()
#include <ctime>   // For seeding random generator


namespace eddycurrent{

Eigen::Matrix<double, 2, 3> GradsBaryCoords(
    Eigen::Matrix<double, 2, 3> vertices) {
  // Compute gradients of barycentric coordinate functions for a flat triangle,
  // whose vertex coordinates are passed in the columns of the argument matrix
  // The algorithm is explained in Remark 2.4.5.9 in the lecture document
  Eigen::Matrix<double, 3, 3> X;
  // See (2.4.5.10) in lecture document: first column of X contains all 1s, the
  // other two columns the first and second coordinates of the vertex coordinate
  // vectors
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // Compute the gradients of the barycentric coordinate functions
  // as columns of a 2x3 matrix containing the \beta-coefficients in (2.4.5.10).
  return X.inverse().block<2, 3>(1, 0);
}

std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<double>, 
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<double>>
  readMeshWithTags(std::string filename, std::map<int, double> tag_to_current, std::map<int, double> tag_to_permeability, std::map<int, double> tag_to_conductivity){
  // std::cout << "Reading mesh from " << filename << std::endl;
  // Total number of contacts
  const int NPhysGrp = 2;
  // load the mesh from a .msh file produced by Gmsh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), filename);
  // Obtain pointer to mesh object
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};
  const lf::mesh::Mesh &mesh{*mesh_p};
  // Output information on the mesh
  //lf::mesh::utils::PrintInfo(std::cout, mesh);
  // A set of integers associated with edges of the mesh (codim = 0 entities)
  lf::mesh::utils::CodimMeshDataSet<double> cell_current{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<double> cell_permeability{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity{mesh_p, 0, -1};
  for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  " edge must be a triangle!");
    cell_current(*cell) = tag_to_current[reader.PhysicalEntityNr(*cell)[0]];
    cell_permeability(*cell) = tag_to_permeability[reader.PhysicalEntityNr(*cell)[0]];
    cell_conductivity(*cell) = tag_to_conductivity[reader.PhysicalEntityNr(*cell)[0]];
  }

  return {mesh_p, cell_current, cell_permeability, cell_conductivity};
}

const Eigen::Matrix<double, 3, 3>  ElemMatProvider::Eval(const lf::mesh::Entity &cell){
  double mu = cell_permeability_(cell);
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  Eigen::Matrix <double , 2 , 3> X = GradsBaryCoords(V);
  Eigen::Matrix <double, 2, 3> temp = X;
  temp.block<1, 3>(0, 0) = X.block<1,3> (1, 0);
  temp.block<1, 3>(1, 0) = X.block<1,3> (0, 0);
  double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  // std::cout << "element matrix\n" << X.transpose() * temp << std::endl;
  return 1/mu * area * temp.transpose() * temp; 
}

const Eigen::Matrix<double, 3, 1> ElemVecProvider::Eval(const lf::mesh::Entity &cell){
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  double cell_current = cell_current_(cell);
  // std::cout << "cell current: " << cell_current << std::endl;
  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  return cell_current * area / 3 * Eigen::Vector3d::Ones();
}

const Eigen::Matrix<double, 3, 3> MassMatProvider::Eval(const lf::mesh::Entity &cell){
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  // std::cout << "cell current: " << cell_current << std::endl;
  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  // std::cout << "element mass matrix" <<  std::endl << cell_conductivity * area / 3 * Eigen::MatrixXd::Identity(3,3) << std::endl;
  // lf::mesh::utils::MeshFunctionConstant gamma(1); 	
  // lf::fe::MassElementMatrixProvider massElemMat(fe_space_, gamma);
  //return massElemMat.Eval(cell) * 2;
  double cell_conductivity = cell_conductivity_(cell);
  return cell_conductivity * area / 3 * Eigen::MatrixXd::Identity(3,3); 
}


Eigen::VectorXd solve(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_permeability, 
    lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity
    ) {

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

  return sol_vec;
}  // end solveMixedBVP

} // namespace eddycurrent