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
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<unsigned>>
  readMeshWithTags(std::string filename, std::map<int, double> tag_to_current, std::map<int, double> tag_to_permeability, std::map<int, double> tag_to_conductivity){
  // Total number of contacts
  const int NPhysGrp = 2;
  // load the mesh from a .msh file produced by Gmsh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), filename);
  // Obtain pointer to mesh object
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};
  const lf::mesh::Mesh &mesh{*mesh_p};
  // Output information on the mesh
  // A set of integers associated with edges of the mesh (codim = 0 entities)
  lf::mesh::utils::CodimMeshDataSet<double> cell_current{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<double> cell_permeability{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tag{mesh_p, 0, 0};
  for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  " edge must be a triangle!");
    cell_current(*cell) = tag_to_current[reader.PhysicalEntityNr(*cell)[0]];
    cell_permeability(*cell) = tag_to_permeability[reader.PhysicalEntityNr(*cell)[0]];
    cell_conductivity(*cell) = tag_to_conductivity[reader.PhysicalEntityNr(*cell)[0]];
    cell_tag(*cell) = reader.PhysicalEntityNr(*cell)[0];
  }
  return {mesh_p, cell_current, cell_permeability, cell_conductivity, cell_tag};
}

const Eigen::Matrix<double, 3, 3>  ElemMatProvider::Eval(const lf::mesh::Entity &cell){
  int material_tag = material_tags_(cell);  
  // Create material if we haven't seen this tag before
  if (materials_.find(material_tag) == materials_.end()) {
      materials_[material_tag] = MaterialFactory::Create(material_tag);
  }
  const auto& material = materials_[material_tag];


  // Get the magnetic flux density and compute reluctivity
  Eigen::Vector2d center_of_triangle;
  center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                   //angle, but I'll take [0.5, 0.5]
    auto magnetic_flux_result = cell_magnetic_flux_(cell, center_of_triangle);
  if (magnetic_flux_result.size() == 0) {
      throw std::runtime_error("cell_magnetic_flux_ returned empty result");
  }
  Eigen::VectorXd B_field = magnetic_flux_result[0];

  double reluctivity = material->getReluctivity(B_field.norm());
  double reluctivity_derivative = material-> getReluctivityDerivative(B_field.norm());



  // std::cout << "reluctivity " << reluctivity << std::endl; 
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  Eigen::Matrix <double , 2 , 3> X = GradsBaryCoords(V);
  Eigen::Matrix <double, 2, 3> temp = X;
  temp.block<1, 3>(0, 0) = X.block<1,3> (1, 0);
  temp.block<1, 3>(1, 0) = X.block<1,3> (0, 0);
  double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  // std::cout << "element matrix A : " << std::endl << reluctivity * area * temp.transpose() * temp << std::endl; 

  // if (material_tag == 3 && B_field.norm() > 0.001){
  //   std::cout << "GradBaryCoords: " << GradsBaryCoords(V) << std::endl;
  //   std::cout << "reluctivity " << reluctivity << std::endl ; 
  //   std::cout << "reluctivity derivative : " << reluctivity_derivative << std::endl ; 
  //   std::cout << "magnetic field " << B_field.norm() << std::endl; 
  //   double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  //   std::cout << "area : " << area << std::endl;
  //   std::cout << "element matrix A : " << std::endl << reluctivity * area * temp.transpose() * temp << std::endl; 
  //   std::cout << std::endl << std::endl; 
  //   // assert(false); 
  // }
  return reluctivity * area * temp.transpose() * temp; 
}



const Eigen::Matrix<double, 3, 3>  ElemMat_N_Provider::Eval(const lf::mesh::Entity &cell){
    int material_tag = material_tags_(cell);  
  // Create material if we haven't seen this tag before
  if (materials_.find(material_tag) == materials_.end()) {
      materials_[material_tag] = MaterialFactory::Create(material_tag);
  }
  
  const auto& material = materials_[material_tag];
  
  Eigen::Vector2d center_of_triangle;
    center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                    //angle, but I'll take [0.5, 0.5
  auto magnetic_flux = cell_magnetic_flux_(cell, center_of_triangle);
  Eigen::VectorXd B_field = magnetic_flux[0];


  double reluctivity_derivative = material->getReluctivityDerivative(B_field.norm());  
  auto grad_xn_result = cell_grad_xn_(cell, center_of_triangle);
  if (grad_xn_result.size() == 0) {
      throw std::runtime_error("cell_grad_xn_ returned empty result");
  }
  Eigen::VectorXd grad_xn = grad_xn_result[0];
  // Check grad_xn dimensions
  if (grad_xn.size() != 2) {  // Assuming it should be 2D
      throw std::runtime_error("Unexpected gradient dimension");
  }
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  Eigen::Matrix <double , 2 , 3> temp = GradsBaryCoords(V);
  Eigen::MatrixXd temp_transpose = temp.transpose() * grad_xn; 

  // if (material_tag == 3 && B_field.norm() > 0.8){
  //   std::cout << "B : " << B_field.norm() << std::endl;
  
  //   std::cout << "grad_xn " << std::endl << grad_xn << std::endl;
  //   std::cout << "GradBaryCoords: " << GradsBaryCoords(V) << std::endl;
  //   std::cout << "reluctivity derivative " << reluctivity_derivative << std::endl ; 
  //   std::cout << "magnetic field " << B_field.norm() << std::endl; 
  //   double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  //   std::cout << "area : " << area << std::endl;
  //   std::cout << "element matrix N : " << std::endl <<  reluctivity_derivative * 2 * area * temp_transpose * temp_transpose.transpose() << std::endl; 
  //   std::cout << std::endl << std::endl; 
  //   // assert(false); 
  // }


  double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));

  return reluctivity_derivative * 2 * area * temp_transpose * temp_transpose.transpose(); 
}






const Eigen::Matrix<double, 3, 1> ElemVecProvider::Eval(const lf::mesh::Entity &cell){
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  double cell_current = cell_current_(cell);
  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  return cell_current * area / 3 * Eigen::Vector3d::Ones();
}

const Eigen::Matrix<double, 3, 3> MassMatProvider::Eval(const lf::mesh::Entity &cell){
  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  double cell_conductivity = cell_conductivity_(cell);
  return cell_conductivity * area / 3 * Eigen::MatrixXd::Identity(3,3); 
}


// Eigen::VectorXd solve(
//     std::shared_ptr<const lf::mesh::Mesh> mesh_p,
//     lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
//     lf::mesh::utils::CodimMeshDataSet<double> & cell_permeability, 
//     lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity
//     ) 
//   {

//   auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
//   const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

//   const lf::base::size_type N_dofs(dofh.NumDofs());
//   LF_ASSERT_MSG(N_dofs == mesh_p -> NumEntities(2),
//                 " N_dofs must agree with number of nodes");

//   auto [A_crs, M_crs, phi] = A_M_phi_assembler(mesh_p, cell_current, cell_permeability, cell_conductivity);

//   const Eigen::SparseMatrix<double> preconditioner_crs = 1 * A_crs + 1e-5 * M_crs;

//   Eigen::VectorXd sol_vec(N_dofs);

//   Eigen::SimplicialLDLT< Eigen::SparseMatrix<double>> preconditioner;
//   preconditioner.compute(preconditioner_crs);

//   sol_vec.setZero();

//   BiCGstab(A_crs, N_dofs, preconditioner, phi, sol_vec, 1e-7, 1000, 2);
//   double rel_res  = 0;
//   rel_res = (A_crs * sol_vec - phi).norm() / phi.norm();
//   std::cout << "relative residuum " << rel_res << std::endl;  

//   return sol_vec;
//  }  // end solveMixedBVP


std::tuple<const Eigen::SparseMatrix<double>, 
           const Eigen::SparseMatrix<double>, 
  Eigen::VectorXd> A_M_phi_assembler
  (
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
  lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags,
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B
  )
{

  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_ASSERT_MSG(N_dofs == mesh_p -> NumEntities(2),
                " N_dofs must agree with number of nodes");

  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);

  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  double alpha = 1;
  ElemMatProvider elemMatProv (cell_tags, cell_B);
  ElemVecProvider elemVecProv (cell_current);
  MassMatProvider massMatProv (cell_conductivity);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMatProv, A);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, massMatProv, M);
  lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);

  auto bd_flags {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

  std::vector<std::pair <long, double>> ess_dof_select {};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
    const lf ::mesh::Entity &dof_node{dofh.Entity(dofnum)}; 
    const Eigen::Vector2d node_pos {
    lf::geometry::Corners(*dof_node.Geometry()).col(0)}; 
    Eigen::Vector2d upper_node;
    upper_node << -5, 0;
    if (bd_flags(dof_node)) {
      ess_dof_select.emplace_back ( true , 0 ) ; 
    } else {  
      ess_dof_select.emplace_back ( false , 0 ) ; 
    }
  }                                 
  Eigen::VectorXd phi_copy = phi;    

  lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &A, &phi](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> {
    return ess_dof_select[dof_idx];}, A, phi);

  lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &M, &phi](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> {
    return ess_dof_select[dof_idx];}, M, phi_copy);

  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  const Eigen::SparseMatrix<double> M_crs = M.makeSparse();
  return {A_crs, M_crs, phi};
}


 std::tuple<const Eigen::SparseMatrix<double>, 
                  Eigen::VectorXd> N_rho_assembler
  (
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags, 
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B,
  lf::fe::MeshFunctionGradFE<double, double>  & cell_grad_xn
  ){
    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs(dofh.NumDofs());
    lf::assemble::COOMatrix<double> N(N_dofs, N_dofs);
    ElemMat_N_Provider elemMat_N_provider(cell_tags, cell_B, cell_grad_xn);
    ElemVec_rho_Provider elemVec_rho_provider(cell_B, cell_tags);
  
    Eigen::VectorXd rho(N_dofs);
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMat_N_provider, N);
    lf::assemble::AssembleVectorLocally(0, dofh, elemVec_rho_provider, rho);

    auto bd_flags {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

    std::vector<std::pair <long, double>> ess_dof_select {};
    for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
      const lf ::mesh::Entity &dof_node{dofh.Entity(dofnum)}; 
      if (bd_flags(dof_node)) {
        ess_dof_select.emplace_back ( true , 0 ) ; 
      } else {  
        ess_dof_select.emplace_back ( false , 0 ) ; 
      }
    }                                 

    // std::cout << "N norm : " << (N.makeSparse()).norm() << std::endl; 
    lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &N, &rho](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> { // probably need to change this 
                                                                                                                                          // since rhs is a combination of different stuff
      return ess_dof_select[dof_idx];}, N, rho);
    // std::cout << "N norm : " << (N.makeSparse()).norm() << std::endl; 

  const Eigen::SparseMatrix<double> N_crs = N.makeSparse();
  return{N_crs, rho}; 
}




const Eigen::Matrix<double, 3, 1> ElemVec_rho_Provider::Eval(const lf::mesh::Entity &cell){

  int material_tag = material_tags_(cell);

  const lf ::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::MatrixXd V = lf::geometry::Corners(*geo_ptr);
  Eigen::Matrix <double , 2 , 3> X = GradsBaryCoords(V);


  if (materials_.find(material_tag) == materials_.end()) {
      materials_[material_tag] = MaterialFactory::Create(material_tag);
  }
  
  const auto& material = materials_[material_tag];
  
  // Get the magnetic flux density and compute reluctivity
  Eigen::Vector2d center_of_triangle;
  center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                   //angle, but I'll take [0.5, 0.5]
  Eigen::VectorXd B_field = cell_magnetic_flux_(cell, center_of_triangle)[0];
  Eigen::Vector2d H_field = material->getH(B_field);
  double H_x = H_field[0];
  double H_y = H_field[1]; 

  Eigen::Vector3d result;

  Eigen::MatrixXd X_transpose = X.transpose();

  result << H_x * X_transpose(0, 1) - H_y * X_transpose(0, 0), 
            H_x * X_transpose(1, 1) - H_y * X_transpose(1, 0), 
            H_x * X_transpose(2, 1) - H_y * X_transpose(2, 0); 
  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  return  area / 3 * result;
}


} // namespace eddycurrent