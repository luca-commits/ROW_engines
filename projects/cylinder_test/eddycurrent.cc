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


//this function is used when we already have read the mesh and just need to recompute the currents
//for the increments
lf::mesh::utils::CodimMeshDataSet<double>
getCellCurrent(std::shared_ptr<const lf::mesh::Mesh> mesh_p, 
               std::map<int, double> tag_to_current, 
               lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tag){
  lf::mesh::utils::CodimMeshDataSet<double> cell_current{mesh_p, 0, -1};

  for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
    cell_current(*cell) = tag_to_current[cell_tag(*cell)];
  }

  return cell_current;
}

std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<unsigned>>
  readMeshWithTags(std::string filename, std::map<int, double> tag_to_current, std::map<int, double> tag_to_conductivity){
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
  lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity{mesh_p, 0, -1};
  lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tag{mesh_p, 0, 0};

  // lf::mesh::utils::CodimMeshDataSet<bool> bd_edge_flags{
  //     lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
    LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  " edge must be a triangle!");
    cell_current(*cell) = tag_to_current[reader.PhysicalEntityNr(*cell)[0]];
    cell_conductivity(*cell) = tag_to_conductivity[reader.PhysicalEntityNr(*cell)[0]];
    cell_tag(*cell) = reader.PhysicalEntityNr(*cell)[0];
  }
  return {mesh_p, cell_current, cell_conductivity, cell_tag};
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


  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  const lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double> &pfe =
      fe_precomp_[ref_el.Id()];
  if (!pfe.isInitialized()) {
    // Accident: cell is of a type not covered by finite element
    // specifications or there is no quadrature rule available for this
    // reference element type
    std::stringstream temp;
    temp << "No local shape function information or no quadrature rule for "
            "reference element type "
         << ref_el;
    throw lf::base::LfException(temp.str());
  }

  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");

  // Physical dimension of the cell
  unsigned world_dim = geo_ptr->DimGlobal();
  // Gram determinant at quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      determinants.size() == pfe.Qr().NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << pfe.Qr().NumPoints());
  // Fetch the transformation matrices for the gradients
  const Eigen::MatrixXd JinvT(
      geo_ptr->JacobianInverseGramian(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      JinvT.cols() == 2 * pfe.Qr().NumPoints(),
      "Mismatch " << JinvT.cols() << " <-> " << 2 * pfe.Qr().NumPoints());
  LF_ASSERT_MSG(JinvT.rows() == world_dim,
                "Mismatch " << JinvT.rows() << " <-> " << world_dim);

  // Element matrix
  Eigen::Matrix3d mat(pfe.NumRefShapeFunctions(), pfe.NumRefShapeFunctions());
  mat.setZero();

  // Loop over quadrature points
  for (unsigned k = 0; k < pfe.Qr().NumPoints(); ++k) {
    const double w = pfe.Qr().Weights()[k] * determinants[k];
    // Transformed gradients
    const auto trf_grad(
        JinvT.block(0, 2 * static_cast<Eigen::Index>(k), world_dim, 2) *
        pfe.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.rows(), 2)
            .transpose());
    // Transformed gradients multiplied with coefficient
    const auto alpha_trf_grad(reluctivity * trf_grad);
    mat += w * (trf_grad.adjoint() * alpha_trf_grad
                );
  }
  // std::cout << std::endl << mat  << std::endl; 
  return mat;
}



const Eigen::Matrix<double, 3, 3>  ElemMat_N_Provider::Eval(const lf::mesh::Entity &cell){

  int material_tag = material_tags_(cell);  

  if (materials_.find(material_tag) == materials_.end()) {
      materials_[material_tag] = MaterialFactory::Create(material_tag);
  }

  const auto& material = materials_[material_tag];
  Eigen::Vector2d center_of_triangle;
  center_of_triangle << 0.5 , 0.5; 
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
  auto temp_copy = temp; 
  Eigen::MatrixXd temp_transpose = temp.transpose() * grad_xn; 
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
  double conductivity = cell_conductivity_(cell);
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Physical dimension of the cell
  const unsigned world_dim = geo_ptr->DimGlobal();

  // Get a quadrature rule of sufficiently high degree on the element
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);
  lf::quad::QuadRuleCache qr_cache;
  const lf::quad::QuadRule qr = qr_cache.Get(cell.RefEl(), 15 * sfl->Degree());
  // Metric factors in quadrature nodes
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  // Fetch values of reaction coefficient gamma at quadrature points:
  // Alocated element matrix
  Eigen::MatrixXd mat(sfl->NumRefShapeFunctions(), sfl->NumRefShapeFunctions());
  mat.setZero();
  // Compute the reference shape functions in the quadrature points
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());
  // Loop over quadrature points to evaluate quadrature formula
  for (unsigned k = 0; k < qr.NumPoints(); ++k) {
    const double w = qr.Weights()[k] * determinants[k];
    mat += w * ((conductivity * rsf.col(k)) * (rsf.col(k).adjoint()));
  }


  return mat;
}

Eigen::VectorXd phi_assembler
  ( std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> & cell_current
  )
{
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_ASSERT_MSG(N_dofs == mesh_p -> NumEntities(2),
                " N_dofs must agree with number of nodes");

  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  
  Eigen::VectorXd phi = Eigen::VectorXd::Constant(N_dofs, 0);
  Eigen::VectorXd phi_boundary = Eigen::VectorXd::Constant(N_dofs, 0);
  
  ElemVecProvider elemVecProv (cell_current);
  lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);

  auto bd_flags_temp {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

  std::vector<std::pair <long, double>> ess_dof_select {};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
    const lf ::mesh::Entity &dof_node{dofh.Entity(dofnum)}; 
    if (bd_flags_temp(dof_node)) {
      ess_dof_select.emplace_back ( true , 0 ) ; 
    } else {  
      ess_dof_select.emplace_back ( false , 0 ) ; 
    }
  }                                 
    lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &A, &phi](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> {
    return ess_dof_select[dof_idx];}, A, phi);
  return phi;
}


std::tuple<Eigen::SparseMatrix<double>, 
           Eigen::SparseMatrix<double>, 
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
  
  Eigen::VectorXd phi = Eigen::VectorXd::Constant(N_dofs, 0);
  Eigen::VectorXd phi_boundary = Eigen::VectorXd::Constant(N_dofs, 0);
  
  ElemMatProvider elemMatProv (cell_tags, cell_B, fe_space);
  ElemVecProvider elemVecProv (cell_current);
  MassMatProvider massMatProv (cell_conductivity, fe_space);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMatProv, A);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, massMatProv, M);
  lf::assemble::AssembleVectorLocally(0, dofh, elemVecProv, phi);


  auto bd_flags_temp {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

  std::vector<std::pair <long, double>> ess_dof_select {};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
    const lf ::mesh::Entity &dof_node{dofh.Entity(dofnum)}; 
    if (bd_flags_temp(dof_node)) {
      ess_dof_select.emplace_back ( true , 0 ) ; 
    } else {  
      ess_dof_select.emplace_back ( false , 0 ) ; 
    }
  }

  lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &A, &phi](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> {
    return ess_dof_select[dof_idx];
  }, A, phi);

  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  const Eigen::SparseMatrix<double> M_crs = M.makeSparse();

  // Print number of non-zero elements

  return {A_crs, M_crs, phi};
}


 std::tuple<Eigen::SparseMatrix<double>, 
                  Eigen::VectorXd>
  N_rho_assembler
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

  const Eigen::SparseMatrix<double> N_crs = N.makeSparse();  
  return{N_crs, rho}; 
}

Eigen::SparseMatrix<double>
  N_assembler
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
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elemMat_N_provider, N);

    auto bd_flags_temp {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

    std::vector<std::pair <long, double>> ess_dof_select {};
    for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
      const lf ::mesh::Entity &dof_node{dofh.Entity(dofnum)}; 
      if (bd_flags_temp(dof_node)) {
        ess_dof_select.emplace_back ( true , 0 ) ; 
      } else {  
        ess_dof_select.emplace_back ( false , 0 ) ; 
      }
    }

    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs); 

    lf::assemble::FixFlaggedSolutionComponents([&ess_dof_select, &N, &phi](lf::assemble::glb_idx_t dof_idx) -> std::pair <bool, double> {
      return ess_dof_select[dof_idx];
    }, N, phi);

    const Eigen::SparseMatrix<double> N_crs = N.makeSparse();  
    return N_crs; 
}


Eigen::VectorXd 
  rho_assembler
(
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags, 
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B
){
    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs(dofh.NumDofs());
    ElemVec_rho_Provider elemVec_rho_provider(cell_B, cell_tags);
    Eigen::VectorXd rho(N_dofs);
    lf::assemble::AssembleVectorLocally(0, dofh, elemVec_rho_provider, rho);

    auto bd_flags_temp {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};

    for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
      rho[dofnum] = 0;
    }

    return rho;
}


// This vector is needed for the Newton iteration
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
  center_of_triangle << 0.5 , 0.5; 

  Eigen::VectorXd B_field = cell_magnetic_flux_(cell, center_of_triangle)[0];
  Eigen::Vector2d H_field = material->getH(B_field);

  double H_x = H_field[0];
  double H_y = H_field[1]; 

  // std::cout << "B field " << B_field.norm() << std::endl;
  // std::cout << "H_x " << H_x << std::endl; 
  // std::cout << "H_y " << H_y << std::endl; 

  Eigen::Vector3d result;

  Eigen::MatrixXd X_transpose = X.transpose();

  result << H_x * X_transpose(0, 1) - H_y * X_transpose(0, 0), 
            H_x * X_transpose(1, 1) - H_y * X_transpose(1, 0), 
            H_x * X_transpose(2, 1) - H_y * X_transpose(2, 0); 

  double area =  0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  return  area / 3 * result;
}

} // namespace eddycurrent