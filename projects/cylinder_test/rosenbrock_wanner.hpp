#ifndef ROSENBROCK_WANNER
#define ROSENBROCK_WANNER

#include <Eigen/Dense> 
#include "increments_no_rotations.hpp"

std::tuple<Eigen::VectorXd, Eigen::VectorXd>  row_step(double timestep, 
                        double t_0, 
                        Eigen::VectorXd current_timestep, 
                        Eigen::SparseMatrix<double> jacobian,
                        Eigen::SparseMatrix<double> M, 
                        Eigen::SparseMatrix<double> M_preconditioner, 
                        Eigen::MatrixXd time_derivative, 
                        std::function<Eigen::VectorXd(double, Eigen::VectorXd, unsigned i)> f,
                        std::shared_ptr<const lf::mesh::Mesh> mesh_p
                        ){

    std::vector<double> b = {2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00 , 4.3586652150845900e-01};
    std::vector<double> b_hat = {3.7810903145819369e-01, -9.6042292212423178e-02, 5.0000000000000000e-01 , 5.0000000000000000e-01};
    std::vector<Eigen::VectorXd> increments = no_rotation::increments(timestep, t_0, current_timestep,
                                                                      jacobian, M, M_preconditioner, time_derivative, f, mesh_p);
    Eigen::VectorXd next_timestep_high = current_timestep;
    Eigen::VectorXd next_timestep_low = current_timestep;
    
    for(unsigned i = 0; i < increments.size(); ++i){
        // std::cout << "increment i norm: " << i << " " << increments[i].norm() << std::endl; 
        next_timestep_high += b[i] * increments[i];
        next_timestep_low += b_hat[i] * increments[i]; 
    }
    return {next_timestep_high, next_timestep_low}; 
    //return {next_timestep_high, next_timestep_low}; 
}

#endif //ROSENBROCK_WANNER
