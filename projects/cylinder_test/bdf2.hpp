#ifndef BDF2_HPP
#define BDF2_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 

#include <lf/fe/fe.h>
#include "BICG_stab.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace bdf2{

    Eigen::VectorXd newton_step(Eigen::SparseMatrix<double>& N,
                                Eigen::SparseMatrix<double>& A, 
                                Eigen::SparseMatrix<double>& M, 
                                double timestep, 
                                Eigen::VectorXd & previous_newton_step, 
                                Eigen::VectorXd & previous_timestep, 
                                Eigen::VectorXd & current_timestep, 
                                Eigen::VectorXd & phi,
                                Eigen::SparseMatrix<double> M_preconditioner
                            ) {
                                    
        lf::assemble::dim_t N_dofs = phi.size();

        Eigen::SparseMatrix<double> lhs = M + 2./3. * timestep * (A + N); 
        Eigen::SparseMatrix<double> preconditioner_matrix =  M_preconditioner + 2./3. * timestep * (A + N); 

        std::cout << "A non zeros" << A.nonZeros() << std::endl;
        std::cout << "M_preconditioner non zeros" << M_preconditioner.nonZeros() << std::endl;
        std::cout << "lhs non zeros" << (preconditioner_matrix).nonZeros() << std::endl;

        Eigen::SparseLU<Eigen::SparseMatrix<double>> preconditioner;
        std::cout << "computing LU decomposition of preconditioner" << std::endl;
        preconditioner.compute(preconditioner_matrix);
        Eigen::VectorXd rhs = (M + 2./3. * timestep * A) * previous_newton_step - 4./3. * M * current_timestep + 1./3. * M * previous_timestep - 2./3. * timestep * phi;
        Eigen::VectorXd next_newton_step(N_dofs);
        next_newton_step.setZero();

        BiCGstab(lhs, N_dofs, preconditioner, rhs, next_newton_step, 1e-7, 100, 1);

        double rel_res = 0.0;

        rel_res = (lhs * next_newton_step - rhs).norm() / rhs.norm();
        std::cout << "relative residuum direct solver: " << rel_res << std::endl;


        if (rhs.norm() > 1e-15) { 
            LF_ASSERT_MSG(rel_res < 1e-4, "Solver failed, residual is greater than 1e-7");
        } else {
            std::cout <<"Solution vector is 0" << std::endl;
        }
        rel_res = (lhs * next_newton_step - rhs).norm() / rhs.norm();

        return next_newton_step;
    }

}
#endif //BDF2_HPP