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
        infile.close();
    }
    else{
        std::cerr << "Unable to open the file!" << std::endl;
    }

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
  
    for (unsigned i = 1; i < timesteps; ++i){

        double rel_angle = angle_step * i;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh",  M_PI * rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
        std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
        unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
        std::cout<< "Final mesh: " << final_mesh << std::endl;
 
        std::cout << "timestep: " << i << std::endl;
        std::string vtk_filename = std::string("vtk_files/time_dependent/rotating/" + mesh_name  + "/eddy_solution_transient_") + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
        // tag_to_current = {{1,0},  {2, 1}, {3, 0}, {4, 0}}; 
        auto [mesh_p, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_permeability, tag_to_conductivity);
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_permeability, cell_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        std::cout << "N dofs: " << dofh.NumDofs() << std::endl;
        Eigen::VectorXd current_timestep_extended = Eigen::VectorXd::Zero(dofh.NumDofs()); 
        std::cout << "current_timestep.norm() " << current_timestep.norm() << std::endl;        
        current_timestep_extended.head(number_stable_dofs) = current_timestep; 
        std::cout << "current_timestep_extended.norm() " << current_timestep_extended.norm() << std::endl;

        Eigen::VectorXd next_timestep = implicit_euler_step(A, M, step_size, current_timestep_extended, phi);

        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        Eigen::VectorXd discrete_solution = next_timestep;

        // discrete_solution[3867] = 1e-7;

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
        
        Eigen::VectorXd backwards_difference = (next_timestep  - current_timestep_extended);
        auto current = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < backwards_difference.rows(); global_idx++) {
            if (node_conductivity->DefinedOn(dofh.Entity(global_idx))){
                current->operator()(dofh.Entity(global_idx)) = backwards_difference[global_idx] * (*node_conductivity)(dofh.Entity(global_idx));
            }
        }

        lf::mesh::utils::CodimMeshDataSet<double> relative_permeability{mesh_p, 0, -1};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            relative_permeability(*cell) = cell_permeability(*cell)/MU_0;
        }

        lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, discrete_solution);
        utils::MeshFunctionCurl2DFE mf_curl(mf_grad);


        utils::MeshFunctionH mf_H(mf_curl, cell_permeability);

        vtk_writer.WritePointData("induced-current", *current);
        std::cout <<"Backwards difference: " <<  backwards_difference.lpNorm<Eigen::Infinity>() << std::endl; 

        vtk_writer.WriteCellData("B", mf_curl);
        vtk_writer.WriteCellData("Conductivity", cell_conductivity);
        vtk_writer.WriteCellData("Prescribed_Current", cell_current);
        vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
        vtk_writer.WriteCellData("H_field", mf_H);

        std::cout << "current_timestep.size() " << current_timestep.size() << std::endl;
        std::cout << "next_timestep.size() " << next_timestep.size() << std::endl;
        current_timestep = next_timestep.head(number_stable_dofs);
    }
    return 0; 
}