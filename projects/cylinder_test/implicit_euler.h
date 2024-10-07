#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H


#include <eigen3>

Eigen::VectorXd implicit_euler_step(Eigen::MatrixXd A, Eigen::MatrixXd B, double timestep, Eigen::VectorXd current_step, Eigen::VectorXd load_vector){
    Eigen::MatrixXd lhs = timestep * A + M;
    Eigen::MatrixXd rhs = M * current_step + timestep * load_vector; 
    
}

#endif //IMPLICIT_EULER_H