#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 

#include <lf/fe/fe.h>
#include "BICG_stab.hpp"
#include <Eigen/IterativeLinearSolvers>


// this function is only used for linear materials, basically an euler step without 
// 
Eigen::VectorXd implicit_euler_step(const Eigen::SparseMatrix<double>& A,  
                                    const Eigen::SparseMatrix<double>& M, 
                                    double timestep, 
                                    const Eigen::VectorXd & current_step,  
                                    const Eigen::VectorXd & load_vector,
                                    const Eigen::SparseMatrix<double> M_preconditioner){

    // Eigen::SparseMatrix<double> lhs = timestep * (A) + M;
    Eigen::SparseMatrix<double> lhs = timestep * (A);


    Eigen::VectorXd rhs = timestep * load_vector;// M * current_step +


    lf::assemble::dim_t N_dofs = load_vector.size();
    Eigen::VectorXd next_timestep(N_dofs);
    next_timestep.setZero();

    // Eigen::SparseMatrix<double> preconditioner_matrix = (A * timestep + M_preconditioner);

    Eigen::SparseMatrix<double> preconditioner_matrix = (A * timestep + M);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> preconditioner;
    std::cout << "computing LU decomposition of preconditioner" << std::endl;
    preconditioner.compute(preconditioner_matrix);

    BiCGstab(lhs, N_dofs, preconditioner, rhs, next_timestep, 1e-7, 1000, 0);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(lhs);
    // next_timestep = solver.solve(rhs);

    double rel_res = 0.0;
    if (rhs.norm() > 1e-17) {
        rel_res = (lhs * next_timestep - rhs).norm() / rhs.norm();
        LF_ASSERT_MSG(rel_res < 1e-4, "Solver failed, residual is greater than 1e-7");
    } else {
        std::cout << " Solution vector is 0 " << std::endl;
    }
    std::cout << "relative residuum euler step: " << rel_res << std::endl; 
    return next_timestep;
}



Eigen::VectorXd newton_step(Eigen::SparseMatrix<double>& N,
                            Eigen::SparseMatrix<double>& A, 
                            Eigen::SparseMatrix<double>& M, 
                            double timestep, 
                            Eigen::VectorXd & previous_newton_step, 
                            Eigen::VectorXd & previous_time_step,
                            Eigen::VectorXd & phi,
                            Eigen::SparseMatrix<double> M_preconditioner
                           ) {
                                
    lf::assemble::dim_t N_dofs = phi.size();


    Eigen::SparseMatrix<double> lhs = (A * timestep + M ) + N * timestep; 
    Eigen::SparseMatrix<double> preconditioner_matrix = (A * timestep + M_preconditioner + N * timestep);
    

    Eigen::SparseLU<Eigen::SparseMatrix<double>> preconditioner;
    std::cout << "computing LU decomposition of preconditioner" << std::endl;
    preconditioner.compute(preconditioner_matrix);

    Eigen::VectorXd rhs = - phi * timestep - M * previous_time_step + M * previous_newton_step + A * previous_newton_step * timestep;

    Eigen::VectorXd next_newton_step(N_dofs);
    next_newton_step.setZero();

    BiCGstab(lhs, N_dofs, preconditioner, rhs, next_newton_step, 1e-7, 100, 1);  

    double rel_res = 0.0;
    rel_res = (lhs * next_newton_step - rhs).norm() / rhs.norm();


    if (rhs.norm() != 0) { 
        
        LF_ASSERT_MSG(rel_res < 1e-4, "Solver failed, residual is greater than 1e-7");
    } else {
        std::cout <<"Solution vector is 0" << std::endl;
    }
    rel_res = (lhs * next_newton_step - rhs).norm() / rhs.norm();

    return next_newton_step;
}

#endif //IMPLICIT_EULER_H