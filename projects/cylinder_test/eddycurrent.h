#ifndef EDDYCURRENT_H_
#define EDDYCURRENT_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

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
  ElemMatProvider(lf::mesh::utils::CodimMeshDataSet<double> cell_permeability) :
                  cell_permeability_(cell_permeability){}
  // Select cells taken into account during cell-oriented assembly
  virtual bool isActive(const lf ::mesh:: Entity & /*cell*/) { return true; }
  // Compute element matrix for a cell, here a fixed-size matrix
  const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity& cell);
  private:
    lf::mesh::utils::CodimMeshDataSet<double> cell_permeability_;
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
  MassMatProvider(lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity) :
                  cell_conductivity_(cell_conductivity){}
    bool isActive(const lf::mesh::Entity& cell){return true;}; 
    const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &cell);
  private:
    lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_;
};



std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<double>, 
          lf::mesh::utils::CodimMeshDataSet<double>,
          lf::mesh::utils::CodimMeshDataSet<unsigned>>
readMeshWithTags(std::string filename, std::map<int, double>, std::map<int, double>, std::map<int, double>);


Eigen::VectorXd solve(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &);

std::tuple<const Eigen::SparseMatrix<double>, 
           const Eigen::SparseMatrix<double>, 
           Eigen::VectorXd
          > A_M_phi_assembler(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &
  );

} //namespace eddycurrent



#endif //EDDYCURRENT_H_