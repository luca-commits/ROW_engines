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
    std::string type_of_rotor = "cylinder_smooth_ring/";
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor;
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
    mergeEverything(mesh_path +  "rotor.msh", mesh_path + "stator.msh", mesh_path + "airgap.msh", final_mesh); 

    std::map<int, double>      tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring
    std::map<int, double> tag_to_permeability{{1,1.00000037 * MU_0}, {2,0.999994 * MU_0}, {3, 0.999994 * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}};

    auto [mesh_p_temp, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_permeability, tag_to_conductivity);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);


    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());

    std::size_t number_stable_dofs = utils::computeDofsWithoutRing(mesh_p_temp, cell_tag, fe_space_temp);

    std::cout << "number_stable_dofs " << number_stable_dofs << std::endl;

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(number_stable_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;


    double angle_step = M_PI / 360 ; 


  
    for (unsigned i = 1; i < timesteps; ++i){

        double rel_angle = angle_step * i;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh",  M_PI * rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
        std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
        std::cout << "Mesh Path :" << final_mesh << std::endl;
        mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
        std::cout<< "Final mesh: " << final_mesh << std::endl;
 
        std::cout << "timestep: " << i << std::endl;
        std::string vtk_filename = std::string("vtk_files/time_dependent/rotating/test_dof_mapping/test_dof_mapping_") + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}}; 
        auto [mesh_p, cell_current, cell_permeability, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_permeability, tag_to_conductivity);
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_permeability, cell_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        std::cout << "N dofs: " << dofh.NumDofs() << std::endl;

        Eigen::VectorXd current_timestep_extended = Eigen::VectorXd::Zero(dofh.NumDofs()); 
        current_timestep_extended[6000] = 1;

        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        Eigen::VectorXd discrete_solution = current_timestep_extended; 

        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < discrete_solution.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = current_timestep_extended[global_idx];
        }
        vtk_writer.WritePointData("test", *nodal_data);

        
    }
    return 0; 
}