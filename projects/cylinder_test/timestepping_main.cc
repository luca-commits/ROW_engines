#include "implicit_euler.h"
#include "eddycurrent_time_dependent.h"
#include "eddycurrent.h"
#include "utils.h"

#include <fstream>

int main (int argc, char *argv[]){

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
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/") + mesh_name;
    std::cout << mesh_path;

    auto time_to_current = [max_current, ramp_up_time](double time){
        double rate = 1/ramp_up_time;
        return time < ramp_up_time ? rate * time * max_current : max_current;
    };

    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    std::cout << "Reading mesh from " << mesh_path << std::endl;
    lf::io::GmshReader reader_temp(std::move(mesh_factory), mesh_path);
    std::shared_ptr<const lf::mesh::Mesh> mesh_p_temp{reader_temp.mesh()};

    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());
    // std::cout << "Solving BVP with " << N_dofs << " FE dofs" << std::endl;


    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(N_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    std::map<int, double>      tag_to_current{};
    std::map<int, double> tag_to_permeability{{1,1.00000037 * MU_0}, {2,0.999994 * MU_0}, {3, 0.999994 * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}};

  
    for (unsigned i = 0; i < timesteps; ++i){

        std::cout << "timestep: " << i << std::endl;
        std::string vtk_filename = std::string("vtk_files/time_dependent/eddy_solution_transient_") + argv[1] + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}}; //1 -> air, 2-> cylinder, 3 -> ring
        auto [mesh_p, cell_current, cell_permeability, cell_conductivity] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_permeability, tag_to_conductivity);
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_permeability, cell_conductivity);
        Eigen::VectorXd next_timestep = implicit_euler_step(A, M, step_size, current_timestep, phi);


        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        Eigen::VectorXd discrete_solution = current_timestep; 

        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

        std::cout << "N_dofs " << dofh_temp.NumDofs() << std::endl;


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
        
        Eigen::VectorXd backwards_difference = (next_timestep  - current_timestep);
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
        current_timestep = next_timestep;
    }
    return 0; 
}