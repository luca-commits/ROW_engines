#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream> 

#include <lf/fe/fe.h>

Eigen::VectorXd implicit_euler_step(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> M, double timestep, Eigen::VectorXd current_step, Eigen::VectorXd load_vector){
    Eigen::SparseMatrix<double> lhs = timestep * A + M ;
    Eigen::VectorXd rhs = M * current_step + timestep * load_vector; 
    // std::cout << "solving qr" << std::endl;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
    qr.compute(lhs);
    Eigen::VectorXd sol_vec = qr.solve(rhs);
    if(sol_vec.norm() != 0){
        double rel_residual = (lhs * sol_vec - rhs).norm() / sol_vec.norm();
        std::cout << "relative residual : " << rel_residual << std::endl; 
        // LF_ASSERT_MSG(rel_residual < 10e-7,
        //             " solver failed, residual is greater than 10e-7");
    }
    // std::cout << "finished solving qr" << std::endl ;
    std::cout << "backwards difference: " << (sol_vec - current_step).norm() << std::endl;


    return sol_vec;
}


#endif //IMPLICIT_EULER_H