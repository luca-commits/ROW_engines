#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"
#include <fstream>

#include <fstream>

int main (int argc, char *argv[]){

    std::cout << "Hiiii " << std::endl;

    double total_time;
    double step_size;
    std::string mesh_name; 
    double max_current = 1000;
    double conductivity_ring = 5.96e7;
    double ramp_up_time = 0.01;
     unsigned newton_steps;

    std::ofstream residual_file;
    residual_file.open("newton_residuals.txt", std::ios::out | std::ios::app); // Append mode in case of restart

    if (!residual_file.is_open()) {
    throw std::runtime_error("Failed to open newton_residuals.txt for writing");
    }

// Write header in case it's a fresh file
residual_file << "# timestep iteration residual time current\n";
residual_file.flush();

    std::ifstream infile("xinput.txt");
    if (infile.is_open()){
        // Read total_time and timesteps
        if (!(infile >> step_size)){
            std::cerr << "Error reading time step size" << std::endl;
        }
        if (!(infile >> total_time)){  
            std::cerr << "Error reading timesteps" << std::endl;
        }
        if (!(infile >> mesh_name)){
            std::cerr << "Error reading mesh name " << std::endl;
        }
        if (!(infile >> max_current)){
            std::cerr << "Error reading mesh name " << std::endl;
        }
        if (!(infile >> conductivity_ring)){
            std::cerr << "Error reading mesh name " << std::endl;
        }
        if (!(infile >> ramp_up_time)){
            std::cerr << "Error reading mesh name " << std::endl;
        }
        if (!(infile >> newton_steps)){
            std::cerr << "Error reading newton steps " << std::endl;
        }
        infile.close();
    }
    else{
        std::cerr << "Unable to open the file!" << std::endl;
    }

    std::cout << "newton_steps  " << newton_steps << std::endl;

    std::cout << std::endl << "timestep: " << step_size << std::endl;

    std::filesystem::path here = __FILE__;
    std::string type_of_rotor = mesh_name;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor + "/";
    std::cout << mesh_path;

    auto time_to_current = [max_current, ramp_up_time](double time){
        double rate = 1/ramp_up_time;
        return time < ramp_up_time ? rate * time * max_current : max_current;
    };

    double rel_angle = 0;
    remeshAirgap(mesh_path + "airgap.geo",  mesh_path + "airgap.msh",  M_PI * rel_angle);
    rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
    std::string final_mesh = mesh_path + "motor_" + std::to_string(0) + ".msh";
    std::cout << "Mesh Path :" << final_mesh << std::endl;
    unsigned numStableNodes = mergeEverything(mesh_path +  "rotor.msh", mesh_path + "stator.msh", mesh_path + "airgap.msh", final_mesh); 

    std::map<int, double>      tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap
    std::map<int, double> tag_to_permeability{{1,1. * MU_0}, {2,1. * MU_0}, {3, 1. * MU_0}, {4,1. * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}, {4, 0}};

    auto [mesh_p_temp, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);


    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());

    std::size_t number_stable_dofs = utils::computeDofsWithoutAirgap(mesh_p_temp, cell_tag, fe_space_temp);

    std::cout << "number_stable_dofs " << number_stable_dofs << std::endl;
    std::cout << "number_stable_nodes " << numStableNodes << std::endl;
    LF_ASSERT_MSG((number_stable_dofs == numStableNodes), "Something went wrong, the number of stable dofs does not correspond to the number of stableNodes");

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(N_dofs); //number_stable_dofs
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    double angle_step = M_PI / 45; 

    std::vector<std::vector<double>> all_timestep_residuals;

    unsigned total_newton_steps = 0; 
    double time = 0; 

  
    for (unsigned i = 0; time <= total_time; ++i, time += step_size){

        std::vector<double> timestep_residuals;

        double rel_angle = 0;//angle_step * i;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh",  M_PI * rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh" , M_PI * rel_angle);
        std::string final_mesh = mesh_path + "motor.msh";
        // std::cout << "final mesh : " << final_mesh << std::endl; 
        unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh" ,final_mesh);
        std::cout << "numStableDofs " << numStableNodes << std::endl ;
 
        std::string vtk_filename = std::string("vtk_files/time_dependent/non-linear/eddy_solution_transient_non_linear") + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 
        std::cout << "timestep: " << i << std::endl;
        std::cout << "time : " << time << std::endl;
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 


        auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        std::cout << "N dofs: " << dofh.NumDofs() << std::endl;

        // Eigen::VectorXd current_timestep_extended = Eigen::VectorXd::Zero(dofh.NumDofs()); 
        // current_timestep_extended.head(number_stable_dofs) = current_timestep;

        Eigen::VectorXd next_newton_step;
        Eigen::VectorXd current_newton_step = current_timestep; //extended          
        
        double newton_tolerance = 1e-7; 
        double newton_residual = 1000; 

        for (unsigned j = 1; newton_residual >  newton_tolerance; ++j){ //j < newton_steps

            ++total_newton_steps; 

            std::cout << "Newton step " << j << std::endl;
            lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_newton_step);
            utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

            auto [N, rho] = eddycurrent::N_rho_assembler(mesh_p, cell_tag, mf_curl, mf_grad);

            // std::cout << "current newton step - current timestep norm : " << (current_newton_step - current_timestep).norm() << std::endl;


            lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_newton_step); //extended
            utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);
            auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);
            next_newton_step = current_newton_step - newton_step(N, A, M, step_size, current_newton_step, current_timestep, phi); //  current_newton_step
            // Eigen::VectorXd next_newton_step = implicit_euler_step(A, M, step_size, current_timestep, phi);

            lf::fe::MeshFunctionGradFE<double, double> mf_grad_next(fe_space, next_newton_step);
            utils::MeshFunctionCurl2DFE mf_curl_next(mf_grad_next);

            auto [A_next, M_next, phi_next] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_next);

            Eigen::SparseMatrix<double> lhs = (step_size * A_next + M);
            Eigen::VectorXd rhs =  M * current_timestep + step_size * phi; //extended
            newton_residual = (lhs * next_newton_step - rhs).norm() / rhs.norm(); 
            std::cout << "newton residual " << newton_residual << std::endl;
                    residual_file << i << " " 
                     << j << " " 
                     << newton_residual << " "
                     << time << " "
                     << time_to_current(time) << "\n";
        residual_file.flush();

            timestep_residuals.push_back(newton_residual);

            current_newton_step = next_newton_step; 
            std::cout << std::endl;
        }

        Eigen::VectorXd next_timestep = current_newton_step;

        lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_timestep); //extended
        utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);

        auto [A_temp, M_temp, phi_temp] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);

        lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, next_timestep);
        utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl);
        auto lhs = (step_size * A + M);
        auto rhs =  M * current_timestep + step_size * phi; //extended
        double rel_residual = (lhs * next_timestep - rhs).norm()  / rhs.norm(); 
        std::cout << "Right hand side : " << rhs.norm() << std::endl; 
        std::cout << "Left hand side : " << (lhs * next_timestep).norm() << std::endl; 
        if (rhs.norm() > 1e-15) std::cout << "Final newton residual " << rel_residual << std::endl; 


        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        Eigen::VectorXd discrete_solution = next_timestep;
        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < discrete_solution.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = discrete_solution[global_idx];
        }
        vtk_writer.WritePointData("A*", *nodal_data);

        // Eigen::VectorXd backwards_difference = (next_timestep  - current_timestep) / step_size; //extended

        // Avoid subtracting similar numbers directly
        Eigen::VectorXd scale = current_timestep.array().abs().max(next_timestep.array().abs());
        Eigen::VectorXd rel_diff = ((next_timestep.array() / scale.array()) - 
                                (current_timestep.array() / scale.array()));
        Eigen::VectorXd backwards_difference = (rel_diff.array() * scale.array()) / step_size;

        std::cout << "backwards_difference : " << (next_timestep  - current_timestep).norm() << std::endl; //extended
        lf::fe::MeshFunctionFE<double, double> mf_backwards_difference(fe_space, backwards_difference); 

        lf::mesh::utils::CodimMeshDataSet<double> induced_current{mesh_p, 0, -1};

        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {

            Eigen::Vector2d center_of_triangle;
            center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                             //angle, but I'll take [0.5, 0.5
            induced_current(*cell) = - mf_backwards_difference(*cell, center_of_triangle)[0] * cell_conductivity(*cell);
        }


        lf::mesh::utils::CodimMeshDataSet<double> relative_permeability{mesh_p, 0, -1};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            int material_tag = cell_tag(*cell); 
            auto material = MaterialFactory::Create(material_tag); 

            Eigen::Vector2d center_of_triangle;
             center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                             //angle, but I'll take [0.5, 0.5
            auto magnetic_flux = mf_curl(*cell, center_of_triangle);
            Eigen::VectorXd B_field = magnetic_flux[0];

            double reluctivity = material->getReluctivity(B_field.lpNorm<Eigen::Infinity>());

            double permeability = 1 / reluctivity; 

            relative_permeability(*cell) = permeability / 1.256e-6;
        }

        Eigen::Vector2d init_value = Eigen::Vector2d::Zero();
        lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> H{mesh_p, 0, init_value};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            int material_tag = cell_tag(*cell); 
            auto material = MaterialFactory::Create(material_tag); 

            Eigen::Vector2d center_of_triangle;
             center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                             //angle, but I'll take [0.5, 0.5
            auto magnetic_flux = mf_curl(*cell, center_of_triangle);
            Eigen::VectorXd B_field = magnetic_flux[0];
            double reluctivity = material->getReluctivity(B_field.lpNorm<Eigen::Infinity>());

            H(*cell) = B_field * reluctivity; 
        }

        vtk_writer.WriteCellData("B", mf_curl);
        vtk_writer.WriteCellData("Conductivity", cell_conductivity);
        vtk_writer.WriteCellData("Prescribed_Current", cell_current);
        vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
        vtk_writer.WriteCellData("H_field", H);
        vtk_writer.WriteCellData("gradients", mf_grad); 
        vtk_writer.WriteCellData("induced_current", induced_current); 
        
        // current_timestep = current_time_step.head(number_stable_dofs);
        // current_timestep = next_timestep.head(number_stable_dofs);
        current_timestep = next_timestep; 
        std::cout << std::endl; 
    }
    residual_file.close();
    std::cout << "total_newton_steps" << total_newton_steps << std::endl ; 

    return 0; 
}