#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 

#include <lf/fe/fe.h>
#include "BICG_stab.hpp"

Eigen::VectorXd implicit_euler_step(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& M, double timestep, const Eigen::VectorXd & current_step,  const Eigen::VectorXd & load_vector){
    Eigen::SparseMatrix<double> lhs = timestep * A + M;
    std::cout << "lhs.norm() " << lhs.norm() << std::endl;
    Eigen::VectorXd rhs = M * current_step + timestep * load_vector; 
    std::cout << "rhs.norm() " << rhs.norm() << std::endl;
    lf::assemble::dim_t N_dofs = load_vector.size();
    Eigen::VectorXd next_timestep(N_dofs);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(lhs);
    next_timestep = solver.solve(rhs);
    std::cout << "next_timestep.norm() " << next_timestep.norm() << std::endl;

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

Eigen::VectorXd newton_step(const Eigen::SparseMatrix<double>& N, 
                            const Eigen::SparseMatrix<double>& A, 
                            const Eigen::SparseMatrix<double>& M, 
                            double timestep, 
                            const Eigen::VectorXd & previous_newton_step, 
                            const Eigen::VectorXd & previous_time_step,  
                            const Eigen::VectorXd & rho, 
                            const Eigen::VectorXd & phi){

    Eigen::SparseMatrix<double> lhs = A * timestep + M + N * timestep ;
    Eigen::VectorXd rhs = - phi * timestep - M * previous_time_step + M * previous_newton_step + A * previous_newton_step * timestep ;
    lf::assemble::dim_t N_dofs = phi.size();
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