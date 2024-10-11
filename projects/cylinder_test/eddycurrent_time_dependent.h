#ifndef EDDYCURRENT_TIME_DEPENDENT_H_
#define EDDYCURRENT_TIME_DEPENDENT_H_

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

namespace eddycurrent_time_dependent{

// class MassMatProvider{
//   public:
//   MassMatProvider(lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity) :
//                   cell_conductivity_(cell_conductivity){}
//     bool isActive(const lf::mesh::Entity& cell){return cell_conductivity_(cell) != 0;}; 
//     const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &cell);
//   private:
//   lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_;
// };

class MassMatProvider{
  public:
  MassMatProvider(std::shared_ptr< const lf::fe::ScalarFESpace<double>> fe_space, 
                                   lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity) :
                  fe_space_(fe_space), cell_conductivity_(cell_conductivity){}
    bool isActive(const lf::mesh::Entity& cell){return cell_conductivity_(cell) != 0;}; 
    const Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &cell);
  private:
    std::shared_ptr< const lf::fe::ScalarFESpace<double>> fe_space_;
    lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_;
};

class ConductivityMeshFunction{
  public:
    ConductivityMeshFunction(lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity) :
    cell_conductivity_(std::move(cell_conductivity)){}
    auto operator()(const lf::mesh::Entity& e,
                  const Eigen::MatrixXd& local) const  {
    std::vector<double> result;
    for (unsigned i = 0; i < 3; ++i){
      result.push_back(cell_conductivity_(e));
    }
    return result; 
  }
  private:
    lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_;
};




std::tuple<std::shared_ptr<const lf::mesh::Mesh>,
           lf::mesh::utils::CodimMeshDataSet<double>, 
           lf::mesh::utils::CodimMeshDataSet<double>,
           lf::mesh::utils::CodimMeshDataSet<double>
          >
readMeshWithTags(std::string filename, std::map<int, double>, std::map<int, double>, std::map<int, double>);


std::tuple<const Eigen::SparseMatrix<double>, 
           const Eigen::SparseMatrix<double>, 
           Eigen::VectorXd
          > A_M_phi_assembler(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &,
    lf::mesh::utils::CodimMeshDataSet<double> &
  );

} //namespace eddycurrent_time_dependent



#endif //EDDYCURRENT_TIME_DEPENDENT_H_