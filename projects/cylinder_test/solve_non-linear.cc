#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"


#include <fstream>

int main (int argc, char *argv[]){

    std::cout << "Hiiii " << std::endl;

    double total_time;
    unsigned timesteps;
    std::string mesh_name; 
    double max_current = 1000;
    double conductivity_ring = 5.96e7;
    double ramp_up_time = 0.01;
     unsigned newton_steps;

    std::ifstream infile("xinput.txt");
    if (infile.is_open()){
        // Read total_time and timesteps
        if (!(infile >> total_time)){
            std::cerr << "Error reading total_time" << std::endl;
        }
        if (!(infile >> timesteps)){  
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

    double step_size = total_time / double(timesteps); 
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
    std::map<int, double> tag_to_permeability{{1,1.00000037 * MU_0}, {2,0.999994 * MU_0}, {3, 0.999994 * MU_0}, {4,1.00000037 * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}, {4, 0}};

    auto [mesh_p_temp, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_permeability, tag_to_conductivity);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);


    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());

    std::size_t number_stable_dofs = utils::computeDofsWithoutAirgap(mesh_p_temp, cell_tag, fe_space_temp);

    std::cout << "number_stable_dofs " << number_stable_dofs << std::endl;
    std::cout << "number_stable_nodes " << numStableNodes << std::endl;
    LF_ASSERT_MSG((number_stable_dofs == numStableNodes), "Something went wrong, the number of stable dofs does not correspond to the number of stableNodes");

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(number_stable_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    double angle_step = M_PI / 360 ; 

    Eigen::VectorXd current_time_step = Eigen::VectorXd::Zero(dofh_temp.NumDofs());

  
    for (unsigned i = 0; i < timesteps; ++i){

        double rel_angle = angle_step * i;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh",  M_PI * rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
        std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
        unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);

 
        std::string vtk_filename = std::string("vtk_files/time_dependent/non-linear/eddy_solution_transient_non_linear") + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 
        std::cout << "time : " << time << std::endl;
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
        // tag_to_current = {{1,0},  {2, 1}, {3, 0}, {4, 0}}; 
        auto [mesh_p, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_permeability, tag_to_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        std::cout << "N dofs: " << dofh.NumDofs() << std::endl;

        Eigen::VectorXd next_newton_step;
        Eigen::VectorXd current_newton_step = current_time_step; 
        
        double newton_tolerance = 1e-7; 
        double newton_residual = 1000; 

        for (unsigned j = 1; newton_residual >  newton_tolerance && j < 10; ++j){ //j < newton_steps
            // if (j > 20) 
            //     throw std::runtime_error("too many newton iterations, no hope for convergence"); 
            std::cout << "timestep: " << i << std::endl;
            std::cout << "Newton step " << j << std::endl;
            lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_newton_step);
            utils::MeshFunctionCurl2DFE mf_curl(mf_grad);
            auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl);
            auto [N, rho] = eddycurrent::N_rho_assembler(mesh_p, cell_tag, mf_curl, mf_grad);

            next_newton_step = current_newton_step - newton_step(N, A, M, step_size, current_newton_step, current_time_step, rho, phi); //  current_newton_step

            auto lhs = (step_size * A + M);
            auto rhs =  M * current_time_step + step_size * phi;

            Eigen::VectorXd residual = lhs * next_newton_step - rhs;
            newton_residual = (lhs * next_newton_step - rhs).norm()/ rhs.norm(); 

            std::cout << "newton residual " << newton_residual << std::endl; 
            current_newton_step = next_newton_step; 
        }


        Eigen::VectorXd next_timestep = next_newton_step;

        lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, next_newton_step);
        utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);

        
        auto [A_temp, M_temp, phi_temp] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);
        auto lhs = (step_size * A_temp + M_temp);
        auto rhs =  M_temp * current_time_step + step_size * phi_temp;
        Eigen::VectorXd residual = lhs * next_newton_step - rhs;
        double rel_residual = (lhs * next_newton_step - rhs).norm() / rhs.norm(); 


        std::cout << "Final newton residual " << rel_residual << std::endl; 

        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        Eigen::VectorXd discrete_solution = next_timestep;
        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < discrete_solution.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = discrete_solution[global_idx];
        }
        vtk_writer.WritePointData("A*", *nodal_data);
        auto node_conductivity = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (const lf ::mesh::Entity *cell : mesh_p -> Entities(0)){
            if(cell_conductivity(*cell) != 0){
                std::span<const lf::mesh::Entity *const> sub_ent_range = cell -> SubEntities(2);
                for (const lf::mesh::Entity *sub_ent : sub_ent_range) {
                    node_conductivity->operator()(*sub_ent) = cell_conductivity(*cell);
                }
            }
        }
        
        Eigen::VectorXd backwards_difference = (next_timestep  - current_time_step);

        auto current = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < backwards_difference.rows(); global_idx++) {
            if (node_conductivity->DefinedOn(dofh.Entity(global_idx))){
                current->operator()(dofh.Entity(global_idx)) = backwards_difference[global_idx] * (*node_conductivity)(dofh.Entity(global_idx));
            }
        }

        lf::mesh::utils::CodimMeshDataSet<double> relative_permeability{mesh_p, 0, -1};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            int material_tag = cell_tag(*cell); 
            auto material = MaterialFactory::Create(material_tag); 

            Eigen::Vector2d center_of_triangle;
             center_of_triangle << 0.5 , 0.5; //theoretically the gradient should be constant on the triangle. So I could take any point on the tri
                                             //angle, but I'll take [0.5, 0.5
            auto magnetic_flux = mf_curl_temp(*cell, center_of_triangle);
            Eigen::VectorXd B_field = magnetic_flux[0];

            double reluctivity = material->getReluctivity(B_field.norm());
            double permeability = 1/reluctivity; 
            relative_permeability(*cell) = permeability / 1.256e-6;
        }

        lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, discrete_solution);
        utils::MeshFunctionCurl2DFE mf_curl(mf_grad);
        utils::MeshFunctionH mf_H(mf_curl, cell_permeability);
        vtk_writer.WritePointData("induced-current", *current);
        vtk_writer.WriteCellData("B", mf_curl);
        vtk_writer.WriteCellData("Conductivity", cell_conductivity);
        vtk_writer.WriteCellData("Prescribed_Current", cell_current);
        vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
        vtk_writer.WriteCellData("H_field", mf_H);
        vtk_writer.WriteCellData("gradients", mf_grad); 
        // current_timestep = current_time_step.head(number_stable_dofs);
        current_time_step = next_timestep; 
    }
    return 0; 
}