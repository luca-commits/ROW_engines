#ifndef INCREMENTS_NO_ROTATIONS_HPP
#define INCREMENTS_NO_ROTATIONS_HPP

#include <vector> 
#include <Eigen/Dense>
#include "eddycurrent.h"


namespace no_rotation{

Eigen::VectorXd function_evaluation(Eigen::VectorXd current_timestep,
                                    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                                    lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags,
                                    double time, 
                                    double max_current, 
                                    double ramp_up_time, 
                                    int x, 
                                    int i,
                                    double timestep,  
                                    bool debug){         
                        

    auto time_to_current = [max_current, ramp_up_time](double time){
        double rate = 1/ramp_up_time;
        return time < ramp_up_time ? rate * time * max_current : max_current;
    };
    std::map<int, double> tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
    lf::mesh::utils::CodimMeshDataSet<double> cell_current = eddycurrent::getCellCurrent(mesh_p, tag_to_current, cell_tags);

    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_timestep);
    utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

    Eigen::VectorXd rho =  eddycurrent::rho_assembler(mesh_p, cell_tags, mf_curl);
    Eigen::VectorXd load = eddycurrent::phi_assembler(mesh_p, cell_current);

    return load - rho; 

}

std::vector<Eigen::VectorXd> increments(double timestep, 
                                        double t_0, 
                                        Eigen::VectorXd current_timestep,
                                        Eigen::SparseMatrix<double> Jacobian, 
                                        Eigen::SparseMatrix<double> M, 
                                        Eigen::MatrixXd time_derivative,
                                        std::function<Eigen::VectorXd(double, Eigen::VectorXd)> f
                                        ){

    unsigned x = t_0 / timestep; 

    double gamma_ii = 4.3586652150845900e-1;
    std::vector<Eigen::VectorXd> increments;

    std::vector<std::vector<double>> gamma = {{-0.87173304301691801}, {-9.0338057013044082e-01, 5.4180672388095326e-02}, {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01}};
    std::vector<std::vector<double>> alpha = {{8.7173304301691801e-01}, {8.4457060015369423e-01, -1.1299064236484185e-01}, {0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00}};

    std::vector<double> gamma_i; 
    unsigned number_stages = 4; 
    for (unsigned i = 0; i < number_stages; ++i){
        double temp = 0; 
        for (unsigned j  = 0; j < i; ++j){
            temp += gamma[i - 1][j];
        }
        gamma_i.push_back(temp); 
    }


    Eigen::SparseMatrix<double> lhs = M - timestep * gamma_ii * Jacobian;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(lhs); 

    for(unsigned i = 0; i < number_stages; ++i){

        double increment_time = t_0; 

        Eigen::VectorXd increments_sum = Eigen::VectorXd::Zero(M.cols());
        Eigen::VectorXd increment_solution = current_timestep; 

        for(unsigned j = 0; j < i; ++j){
            increments_sum += gamma[i - 1][j] * increments[j];
            increment_time += timestep * alpha[i - 1][j]; 
            increment_solution += increments[j] * alpha[i - 1][j]; 
        }

        Eigen::VectorXd rhs =   timestep * f(increment_time, increment_solution) 
                              + timestep * Jacobian * increments_sum 
                              + timestep * timestep * gamma_i[i] *  time_derivative;



        Eigen::VectorXd increment = solver.solve(rhs); 
        double rel_res = 0.0;
        if (rhs.norm() != 0) {
            rel_res = (lhs * increment - rhs).norm() / rhs.norm();
            LF_ASSERT_MSG(rel_res < 1e-7, "Solver failed, residual is greater than 1e-7");
        } else {
            std::cout << " Solution vector is 0 " << std::endl;
        }
        // std::cout << "relative residuum for incerement " << i << " : " << rel_res << std::endl; 

        increments.push_back(increment);

    }

    return increments;
}
} //namespace no_rotation

#endif //INCREMENTS_NO_ROTATIONS_HPP