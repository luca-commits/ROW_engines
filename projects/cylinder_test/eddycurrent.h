#ifndef EDDYCURRENT_H_
#define EDDYCURRENT_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include "utils.h"
#include "magnetic_material.hpp"

#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include <tuple> 

#define MU_0 4 * M_PI * 1e-7

namespace eddycurrent{

/** @brief Computation of gradients of barycentric coordinate functions. Taken from the homeworks repository
 *         https://github.com/erickschulz/NPDECODES/blob/master/homeworks/StationaryCurrents/mastersolution/stationarycurrents.h
 *
 * @param vertices 2x3 matrix whose columns contain the vertex coordinates of
 * the triangle
 * @return gradients of barycentric coordinate functions stored in the columns
 * of a 2x3 matrix.
 *
 */
Eigen::Matrix<double, 2, 3> GradsBaryCoords(
    Eigen::Matrix<double, 2, 3> vertices);

class ElemMatProvider {
 public:
  // Constructor can be used to pass data required for local computations
  ElemMatProvider(lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags,
                  utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux, 
                  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space) :
                   material_tags_(material_tags), 
                   fe_precomp_(),
                   cell_magnetic_flux_(cell_magnetic_flux){
      for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
        auto fe = fe_space->ShapeFunctionLayout(ref_el);
        // Check whether shape functions for that entity type are available.
        // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
        // object is not initialized if the associated description of local shape
        // functions is missing.
        if (fe != nullptr) {
          // Precompute cell-independent quantities based on quadrature rules
          // with twice the degree of exactness compared to the degree of the
          // finite element space.
          fe_precomp_[ref_el.Id()] = lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
              fe, lf::quad::make_QuadRule(ref_el, 2 * fe->Degree()));
        }
      }
    }
  // Select cells taken into account during cell-oriented assembly
  virtual bool isActive(const lf ::mesh:: Entity & /*cell*/) { return true; }
  // Compute element matrix for a cell, here a fixed-size matrix
  const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity& cell);
  private:
    lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags_;
    utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux_;
    std::map<int, std::shared_ptr<MagneticMaterial>> materials_;
    std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>, 5> fe_precomp_;
};

class ElemVecProvider{
  public:
  ElemVecProvider(lf::mesh::utils::CodimMeshDataSet<double> cell_current) :
                  cell_current_(cell_current){}
    bool isActive(const lf::mesh::Entity& cell){return cell_current_(cell) != 0;}; //maybe change it to only do computations if J != 0 
    const Eigen::Matrix<double, 3, 1> Eval(const lf::mesh::Entity &cell);
  private:
    lf::mesh::utils::CodimMeshDataSet<double> cell_current_;
};


class MassMatProvider{
  public:
  MassMatProvider(lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity, std::shared_ptr<const lf::fe::ScalarFESpace<double>> fe_space) :
                  cell_conductivity_(cell_conductivity), fe_space_(fe_space){}
    bool isActive(const lf::mesh::Entity& cell){return true;}; 
    const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &cell);
  private:
    lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_;
    std::shared_ptr<const lf::fe::ScalarFESpace<double>> fe_space_;
};

class ElemMat_N_Provider {
public:
    ElemMat_N_Provider(
        lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags,
        utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux,  
        lf::fe::MeshFunctionGradFE<double, double> cell_grad_xn)
        : material_tags_(material_tags)
        , cell_grad_xn_(cell_grad_xn)
        , cell_magnetic_flux_(cell_magnetic_flux) {}
        
    virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  // Compute element matrix for a cell, here a fixed-size matrix
  const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity& cell);
private:
    lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags_;
    utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux_;
    lf::fe::MeshFunctionGradFE<double, double> cell_grad_xn_;
    std::map<int, std::shared_ptr<MagneticMaterial>> materials_;
};




class ElemVec_rho_Provider{
  public:
    ElemVec_rho_Provider(utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux,
                  lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags) :
                  cell_magnetic_flux_(cell_magnetic_flux),
                  material_tags_(material_tags) {}
    bool isActive(const lf::mesh::Entity& cell){return true;}; //maybe change it to only do computations if J != 0 
    const Eigen::Matrix<double, 3, 1> Eval(const lf::mesh::Entity &cell);
  private:
    utils::MeshFunctionCurl2DFE<double, double> cell_magnetic_flux_;
    std::map<int, std::shared_ptr<MagneticMaterial>> materials_;
    lf::mesh::utils::CodimMeshDataSet<unsigned int> material_tags_;
};





std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<double>, 
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<unsigned>>
readMeshWithTags(std::string filename, std::map<int, double>, std::map<int, double>);


Eigen::VectorXd solve(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &);

std::tuple<Eigen::SparseMatrix<double>, 
           Eigen::SparseMatrix<double>, 
           Eigen::VectorXd
          > A_M_phi_assembler(
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<double> & cell_current,
  lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags,
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B
  );



   std::tuple<Eigen::SparseMatrix<double>, 
                  Eigen::VectorXd> N_rho_assembler
  (
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags, 
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B,
  lf::fe::MeshFunctionGradFE<double, double>  & cell_grad_xn
  );

} //namespace eddycurrent



#endif //EDDYCURRENT_H_  