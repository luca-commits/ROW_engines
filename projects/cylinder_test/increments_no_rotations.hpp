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
                    double ramp_up_time){                
    auto time_to_current = [max_current, ramp_up_time](double time){
        double rate = 1/ramp_up_time;
        return time < ramp_up_time ? rate * time * max_current : max_current;
    };
    std::map<int, double> tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
    lf::mesh::utils::CodimMeshDataSet<double> cell_current = eddycurrent::getCellCurrent(mesh_p, 
                                                                                        tag_to_current,
                                                                                        cell_tags);

    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    std::cout<< " current timestep norm " << current_timestep.norm() << std::endl; 

    lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_timestep);
    utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

    Eigen::VectorXd rho =  eddycurrent::rho_assembler(mesh_p, cell_tags, mf_curl);
    Eigen::VectorXd load = eddycurrent::phi_assembler(mesh_p, cell_current);

    std::cout <<  "rho.norm() " << rho.norm() << std::endl; 
    std::cout << "load.norm() " << load.norm() << std::endl; 

    return load - rho; 

}

std::vector<Eigen::VectorXd> increments(double timestep, 
                                        double t_0, 
                                        Eigen::VectorXd current_timestep,
                                        Eigen::SparseMatrix<double> Jacobian, 
                                        Eigen::SparseMatrix<double> M, 
                                        Eigen::MatrixXd time_derivative,
                                        std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                                        lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags,
                                        double max_current, 
                                        double ramp_up_time
                                         ){

    double gamma_ii = 4.3586652150845900e-01;
    std::vector<Eigen::VectorXd> increments;
    std::vector<std::vector<double>> gamma = {{-0.87173304301691801}, {-9.0338057013044082e-01, 5.4180672388095326e-02}, {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01}};
    std::vector<std::vector<double>> alpha = {{8.7173304301691801e-01}, {8.4457060015369423e-01, -1.1299064236484185e-01}, {0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00}};
    unsigned number_stages = 4; 

    Eigen::MatrixXd k_0 = Eigen::VectorXd::Zero(M.cols());  

    Eigen::SparseMatrix<double> lhs = M - timestep * gamma_ii * Jacobian;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.compute(lhs); 

    for(unsigned i = 0; i < number_stages; ++i){

        Eigen::VectorXd temp = Eigen::VectorXd::Zero(M.cols());
        double increment_time = 0; 
        Eigen::VectorXd increment_solution = current_timestep; 

        for(unsigned j = 0; j < i; ++j){
            std::cout << "i " << i << " j " << j << std::endl;
            temp += timestep * gamma[i - 1][j] * Jacobian * increments[j];
            increment_time += timestep * alpha[i - 1][j]; 
            increment_solution += increments[j] * alpha[i - 1][j]; 
        }

        Eigen::VectorXd rhs = timestep * timestep * time_derivative + temp + timestep * function_evaluation(increment_solution, mesh_p, cell_tags,  increment_time, max_current, ramp_up_time);
        std::cout << "temp.norm(): " << temp.norm() << std::endl; 
        std::cout << "time derivative norm " << time_derivative.norm() << std::endl; 
        std::cout << "function evaluation norm " << function_evaluation(increment_solution, mesh_p, cell_tags,  increment_time, max_current, ramp_up_time).norm(); 
        std::cout << "rhs.norm" << rhs.norm() << std::endl; 
        std::cout << "increment i norm " << i << " " << solver.solve(rhs).norm() << std::endl; 

        Eigen::VectorXd increment = solver.solve(rhs); 
            double rel_res = 0.0;
        if (rhs.norm() != 0) {
            rel_res = (lhs * next_timestep - rhs).norm() / rhs.norm();
            LF_ASSERT_MSG(rel_res < 1e-4, "Solver failed, residual is greater than 1e-7");
        } else {
            std::cout << " Solution vector is 0 " << std::endl;
        }

        increments.push_back(increment);


        std::string vtk_filename = std::string("vtk_files/time_dependent/debug") + std::to_string(i) + std::string(".vtk");
        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);

        std::cout << "vtk file " << vtk_filename << std::endl; 
        vtk_writer.setBinary(true);


        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < increment.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = increment[global_idx];
        }
        vtk_writer.WritePointData("A*", *nodal_data);

    }

    return increments;
}
} //namespace no_rotation

#endif //INCREMENTS_NO_ROTATIONS_HPP