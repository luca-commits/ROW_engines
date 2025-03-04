#include "rosenbrock_wanner.hpp"
#include <fstream>
#include <Eigen/Sparse>

int main(int argc, char *argv[]) {
    double total_time;
    double step_size;

    std::ifstream infile("xinput_test.txt");
    if (infile.is_open()) {
        if (!(infile >> step_size)) {
            std::cerr << "Error reading time step size" << std::endl;
        }
        if (!(infile >> total_time)) {
            std::cerr << "Error reading timesteps" << std::endl;
        }
    } else {
        std::cerr << "Unable to open the file!" << std::endl;
        return 1;
    }
    infile.close();

    std::ofstream outfile("results.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error opening results file!" << std::endl;
        return 1;
    }
    outfile << "Time,NumSol_0,NumSol_1,ExactSol_0,ExactSol_1\n";

    unsigned N_dofs = 2;
    Eigen::VectorXd initial_timestep = Eigen::VectorXd::Ones(N_dofs);
    Eigen::VectorXd current_timestep = initial_timestep;
    double current_time = 0;

    Eigen::SparseMatrix<double> M(N_dofs, N_dofs);
    M.insert(0, 0) = 1.0;
    M.insert(1, 1) = 1.0;

    Eigen::SparseMatrix<double> jacobian(N_dofs, N_dofs);
    jacobian.insert(0, 0) = -1.0;
    jacobian.insert(1, 1) = -5.0;

    Eigen::VectorXd time_derivative = Eigen::VectorXd::Zero(N_dofs);

    auto function_evaluator = [](double c_time, Eigen::VectorXd current_timestep) -> Eigen::VectorXd {
        Eigen::VectorXd b(2);
        b(0) = -1.0 * current_timestep(0);
        b(1) = -5.0 * current_timestep(1);
        return b;
    };

    for (unsigned i = 1; current_time <= total_time; ++i, current_time += step_size) {
        current_time = i * step_size;
        Eigen::VectorXd next_timestep = row_step(step_size, current_time, current_timestep, jacobian, M, time_derivative, function_evaluator);

        Eigen::VectorXd exact_solution(2);
        exact_solution << initial_timestep[0] * std::exp(-1 * current_time),
                          initial_timestep[1] * std::exp(-5 * current_time);

        double l2_error_norm = (next_timestep - exact_solution).norm();
        std::cout << "l2_error_norm : " << l2_error_norm << std::endl;

        outfile << current_time << "," << next_timestep(0) << "," << next_timestep(1) << ","
                << exact_solution(0) << "," << exact_solution(1) << "\n";

        current_timestep = next_timestep;
    }

    outfile.close();
    return 0;
}
