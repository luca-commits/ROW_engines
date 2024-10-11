#include "eddycurrent.h"
#include "eddycurrent_time_dependent.h"


namespace eddycurrent_time_dependent{

// const Eigen::Matrix<double, 3, 3> MassMatProvider::Eval(const lf::mesh::Entity &cell){
//   const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
//   Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
//   double cell_conductivity = cell_conductivity_(cell);
//   // std::cout << "cell current: " << cell_current << std::endl;
//   double area =  std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
//   // std::cout << "element mass matrix" <<  std::endl << cell_conductivity * area / 3 * Eigen::MatrixXd::Identity(3,3) << std::endl;
//   return cell_conductivity * (area / 3) * Eigen::MatrixXd::Identity(3,3);
// }

const Eigen::Matrix<double, 3, 3> MassMatProvider::Eval(const lf::mesh::Entity &cell){
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  // std::cout << "cell current: " << cell_current << std::endl;
  double area =  std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  // std::cout << "element mass matrix" <<  std::endl << cell_conductivity * area / 3 * Eigen::MatrixXd::Identity(3,3) << std::endl;
  ConductivityMeshFunction MF(cell_conductivity_);
  lf::fe::MassElementMatrixProvider massElemMat(fe_space_, MF);
  return massElemMat.Eval(cell);
}


std::tuple<const Eigen::SparseMatrix<double>, 
           const Eigen::SparseMatrix<double>, 
           Eigen::VectorXd> A_M_phi_assembler
    (
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_permeability, 
    lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity 
    )
{

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
  // std::cout << "phi\n" << phi << std::endl;
  eddycurrent::ElemMatProvider elemMatProv(cell_permeability);
  eddycurrent::ElemVecProvider elemVecProv(cell_current);
  MassMatProvider massMatProv (fe_space, cell_conductivity);

  // Cell-oriented assembly
  // std::cout << "assebling matrix" << std::endl;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMatProv, A);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, massMatProv, M);
  lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);


  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  const Eigen::SparseMatrix<double> M_crs = M.makeSparse();

  // std::cout << "M non zeors " <<  M_crs.nonZeros() << std::endl;

  // std::cout << "M " << std::endl << M_crs << std::endl;

  // std::cout << "COO M " << std::endl << M_crs << std::endl;
  return {A_crs, M_crs, phi};
}

std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<double>, 
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<double>>
readMeshWithTags(
  std::string filename, 
  std::map<int, double> tag_to_current, 
  std::map<int, double> tag_to_permeability, 
  std::map<int, double> tag_to_conductivity){
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
    cell_current(*cell) =      tag_to_current[reader.PhysicalEntityNr(*cell)[0]];
    cell_permeability(*cell) = tag_to_permeability[reader.PhysicalEntityNr(*cell)[0]];
    cell_conductivity(*cell) = tag_to_conductivity[reader.PhysicalEntityNr(*cell)[0]];
  }

  return {mesh_p, cell_current, cell_permeability, cell_conductivity};
}

} // namespace eddycurrent_time_dependent