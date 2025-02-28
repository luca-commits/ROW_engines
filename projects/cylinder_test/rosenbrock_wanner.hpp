#ifndef ROSENBROCK_WANNER
#define ROSENBROCK_WANNER

#include <Eigen/Dense> 
#include "increments_no_rotations.hpp"

Eigen::VectorXd row_step(double timestep, 
                        double t_0, 
                        Eigen::VectorXd current_timestep, 
                        Eigen::SparseMatrix<double> jacobian,
                        Eigen::SparseMatrix<double> M, 
                        Eigen::MatrixXd time_derivative, 
                        std::function<Eigen::VectorXd(double, Eigen::VectorXd)> f
                        ){

    std::vector<double> b = {2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00 , 4.3586652150845900e-01};
    std::vector<Eigen::VectorXd> increments = no_rotation::increments(timestep, t_0, current_timestep,
                                                                      jacobian, M, time_derivative, f);
    Eigen::VectorXd next_timestep = current_timestep;
    for(unsigned i = 0; i < increments.size(); ++i){
        std::cout << "increment i norm: " << i << " " << increments[i].norm() << std::endl; 
        next_timestep += b[i] * increments[i];
    }
    return next_timestep; 
}

#endif //ROSENBROCK_WANNER
