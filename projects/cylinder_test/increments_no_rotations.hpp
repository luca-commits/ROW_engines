#ifndef INCREMENTS_NO_ROTATIONS_HPP
#define INCREMENTS_NO_ROTATIONS_HPP

#include <Eigen/Dense>
#include <vector>

#include "eddycurrent.h"

namespace no_rotation {

std::vector<Eigen::VectorXd> increments(
    double timestep, double t_0, Eigen::VectorXd current_timestep,
    Eigen::SparseMatrix<double> Jacobian, Eigen::SparseMatrix<double> M,
    Eigen::SparseMatrix<double> M_preconditioner,
    Eigen::MatrixXd time_derivative,
    std::function<Eigen::VectorXd(double, Eigen::VectorXd, unsigned i)> f,
    std::shared_ptr<const lf::mesh::Mesh> mesh_p
    ) {
  unsigned x = t_0 / timestep;

  double gamma_ii = 4.3586652150845900e-1; //this is the gamma which defines the method
  std::vector<Eigen::VectorXd> increments;

  std::vector<std::vector<double>> gamma = {
      {-0.87173304301691801},
      {-9.0338057013044082e-01, 5.4180672388095326e-02},
      {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01}};
  std::vector<std::vector<double>> alpha = {
      {8.7173304301691801e-01},
      {8.4457060015369423e-01, -1.1299064236484185e-01},
      {0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00}};

  std::vector<double> gamma_i;
  unsigned number_stages = 4;
  for (unsigned i = 0; i < number_stages; ++i) {
    double temp = 0;
    for (unsigned j = 0; j < i; ++j) {
      temp += gamma[i - 1][j];
    }
    gamma_i.push_back(temp);
  }

  Eigen::SparseMatrix<double> lhs = -timestep * gamma_ii * Jacobian + M;//M - timestep * gamma_ii * Jacobian;
  Eigen::SparseMatrix<double> preconditioner_matrix = -timestep * gamma_ii * Jacobian + M;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> preconditioner;
  std::cout << "computing LU decomposition of preconditioner" << std::endl;
  preconditioner.compute(preconditioner_matrix);

  for (unsigned i = 0; i < number_stages; ++i) {
    double increment_time = t_0;

    Eigen::VectorXd increments_sum = Eigen::VectorXd::Zero(M.cols());
    Eigen::VectorXd increment_solution = current_timestep;

    for (unsigned j = 0; j < i; ++j) {
      increments_sum += gamma[i - 1][j] * increments[j];
      increment_time += timestep * alpha[i - 1][j];
      increment_solution += increments[j] * alpha[i - 1][j];
    }

    Eigen::VectorXd rhs = timestep *  f(increment_time, increment_solution, i) +
                          timestep * Jacobian * increments_sum +
                          timestep * timestep * gamma_i[i] * time_derivative;

    lf::assemble::dim_t N_dofs = rhs.size();
    Eigen::VectorXd increment(N_dofs);
    increment.setZero();
    BiCGstab(lhs, N_dofs, preconditioner, rhs, increment, 1e-16, 100, 1);
    double rel_res = 0.0;
    if (rhs.norm() != 0) {
      rel_res = (lhs * increment - rhs).norm() / rhs.norm();
      LF_ASSERT_MSG(rel_res < 1e-7,
                    "Solver failed, residual is greater than 1e-7");
    } else {
      std::cout << " Solution vector is 0 " << std::endl;
    }
    std::cout << "relative residuum for incerement " << i << " : " << rel_res
              << std::endl;

    increments.push_back(increment);

    bool debug = 0; 
    if (debug){
        std::string vtk_filename = std::string("vtk_files/time_dependent/debug_") + std::to_string(4*x + i) + std::string(".vtk");
        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);

        std::cout << "vtk file " << vtk_filename << std::endl; 
        vtk_writer.setBinary(true);

        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = increment[global_idx];
        }
        vtk_writer.WritePointData("increment", *nodal_data);

        Eigen::VectorXd timestep_jacobian = timestep * Jacobian * increments_sum;

        nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = timestep_jacobian[global_idx];
        }
        vtk_writer.WritePointData("jacobian-increment_sum", *nodal_data);

        Eigen::VectorXd time_derivative_timestep = gamma_i[i]  * timestep * timestep * time_derivative;
        nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx))= time_derivative_timestep[global_idx];
        }
        vtk_writer.WritePointData("timestep_derivative", *nodal_data);

        nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx))= rhs[global_idx];
        }
        vtk_writer.WritePointData("rhs", *nodal_data);

        Eigen::VectorXd feval = timestep * f(increment_time, increment_solution, i);
        nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx))= feval[global_idx];
        }
        vtk_writer.WritePointData("feval", *nodal_data);

        nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx))= increment_solution[global_idx];
        }
        vtk_writer.WritePointData("increment_solution", *nodal_data);
    }

    
  }

  return increments;
}
}  // namespace no_rotation

#endif  // INCREMENTS_NO_ROTATIONS_HPP