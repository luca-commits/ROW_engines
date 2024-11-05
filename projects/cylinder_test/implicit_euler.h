#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 

#include <lf/fe/fe.h>
#include "BICG_stab.hpp"

Eigen::VectorXd implicit_euler_step(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& M, double timestep, const Eigen::VectorXd & current_step, const Eigen::VectorXd & load_vector){
    // Eigen::SparseMatrix<double> Id(load_vector.size(), load_vector.size());
    // Id.setZero(); 
    Eigen::SparseMatrix<double> lhs = timestep * A + M;
    Eigen::VectorXd rhs = M * current_step + timestep * load_vector; 
    lf::assemble::dim_t N_dofs = load_vector.size();
    Eigen::VectorXd next_timestep(N_dofs);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(lhs);
    next_timestep = solver.solve(rhs);


    double rel_res = 0.0;
    if (rhs.norm() != 0) {
        rel_res = (lhs * next_timestep - rhs).norm() / rhs.norm();
        if (rel_res > 1e-7) {
            LF_ASSERT_MSG(rel_res < 1e-7, "Solver failed, residual is greater than 1e-7");
        }
    } else {
        std::cout <<"Solution vector is 0" << std::endl;
    }

    return next_timestep;
}
#endif //IMPLICIT_EULER_H