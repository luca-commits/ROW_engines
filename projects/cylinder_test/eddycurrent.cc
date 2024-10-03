#include "eddycurrent.h"

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
          lf::mesh::utils::CodimMeshDataSet<double>>
  readMeshWithTags(std::string filename, std::map<int, double> tag_to_current, std::map<int, double> tag_to_permeability){
  std::cout << "Reading mesh from " << filename << std::endl;
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
  for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  " edge must be a triangle!");
    cell_current(*cell) = tag_to_current[reader.PhysicalEntityNr(*cell)[0]];
    cell_permeability(*cell) = tag_to_permeability[reader.PhysicalEntityNr(*cell)[0]];
  }

  return {mesh_p, cell_current, cell_permeability};
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
  double area =  std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  return cell_current * area / 3 * Eigen::Vector3d::Ones();
}

Eigen::VectorXd solve(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
     lf::mesh::utils::CodimMeshDataSet<double> & cell_permeability
    ) {
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
  // Right-hand side vector
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  // std::cout << "phi\n" << phi << std::endl;
  ElemMatProvider elemMatProv (cell_permeability);
  ElemVecProvider elemVecProv (cell_current) ;
  // Cell-oriented assembly
  std::cout << "assebling matrix" << std::endl;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMatProv, A);
  lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);

  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  std::cout << "solving qr" << std::endl;
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
  qr.compute(A_crs);
  Eigen::VectorXd sol_vec = qr.solve(phi);
  std::cout << "finished solving qr" << std::endl ;


  // LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // // std::cout << "phi \n" << phi << std::endl;
  // // Eigen::VectorXd sol_vec = solver.solve(phi);
  // LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  //std::cout << "sol vec \n" << sol_vec << std::endl;
  return sol_vec;
}  // end solveMixedBVP

} // namespace eddycurrent