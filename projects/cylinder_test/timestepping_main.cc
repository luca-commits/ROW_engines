#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"

#include <fstream>

int main (int argc, char *argv[]){

    double total_time;
    double step_size;
    std::string mesh_name; 
    double max_current = 1000;
    double conductivity_ring = 5.96e7;
    double ramp_up_time = 0.01;
     unsigned newton_steps;

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

    std::cout << std::endl << "timestep: " << step_size << std::endl;

    std::filesystem::path here = __FILE__;
    mesh_name = "rotating/box_irregular_ring/motor.msh";
    std::cout << "mesh name" <<mesh_name << std::endl;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/") + mesh_name;

    auto time_to_current = [max_current, ramp_up_time](double time){
        double rate = 1/ramp_up_time;
        return time < ramp_up_time ? rate * time * max_current : max_current;
    };


    std::map<int, double> tag_to_conductivity{{1, 1e-10}, {2, 1e-10}, {3, conductivity_ring}, {4, 1e-10}}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap
    std::map<int, double>  tag_to_current = {{1, 0},  {2, 1e-30}, {3, 0}, {4, 0}}; //1 -> air, 2-> cylinder, 3 -> ring

    auto [mesh_p_temp, cell_current_temp, cell_conductivity_temp, cell_tag_temp] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_conductivity, 0);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());
    

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Constant(N_dofs, 0e-11); 


    // Eigen::VectorXd previous_timestep = Eigen::VectorXd::Zero(N_dofs);
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    // lf::mesh::utils::CodimMeshDataSet<bool> current_indicator{mesh_p_temp, 0, false};
    // for (const lf::mesh::Entity *cell : mesh_p_temp->Entities(0)) {
    //     if (cell_current_temp(*cell) > 0) {
    //         current_indicator(*cell) = true;
    //     }
    // }
    // double cylinder_area = utils::computeArea(mesh_p_temp, current_indicator);

    // //create a mesh dataset that is one for cells on boundary
    // lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
    //     lf::mesh::utils::flagEntitiesOnBoundary(mesh_p_temp, 2)};

    // lf::mesh::utils::CodimMeshDataSet<bool> bd_indicator{mesh_p_temp, 0, false};
    // for (const lf::mesh::Entity *cell : mesh_p_temp->Entities(0)) {
    //     bool on_boundary = false;
    //     for (const lf::mesh::Entity *node : cell->SubEntities(2)) {
    //         if (bd_flags(*node)) {
    //             on_boundary = true;
    //             break;
    //         }
    //     }
    //     if (on_boundary) {
    //         bd_indicator(*cell) = true;
    //     }
    // }

    // double boundary_area = utils::computeArea(mesh_p_temp, bd_indicator);

    // std::cout << "boundary_area " << boundary_area << std::endl; 

    // std::cout << "cylinder_area " << cylinder_area << std::endl;

    double time = 0; 

    for (unsigned i = 0; time <= total_time; ++i, time += step_size){


        double boundary_current = time_to_current(time); //ßcylinder_area  / boundary_area * time_to_current(time);



        std::cout << "timestep: " << i << std::endl;
        std::string vtk_filename = std::string("vtk_files/time_dependent/eddy_solution_transient_") + argv[1] + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 

        std::cout << "current : " << time_to_current(time) << std::endl ; 

        std::cout << "time_to_current " << time_to_current(time) << std::endl;
        // std::cout << "boundary_current " << boundary_current << std::endl;

        std::cout << "time: " << time << std::endl; 

        tag_to_current = {{1, 0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; //1 -> air, 2-> cylinder, 3 -> ring
        

        auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(mesh_path, tag_to_current, tag_to_conductivity, boundary_current);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
        // Obtain local->global index mapping for current finite element space
        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        // Dimension of finite element space = number of nodes of the mesh

        lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_timestep);
        utils::MeshFunctionCurl2DFE mf_curl(mf_grad);


        auto [A, M, phi, phi_boundary] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl);  
        std::cout << "phi boundary norm: " << phi_boundary.norm() << std::endl; 
        Eigen::VectorXd next_timestep = implicit_euler_step(A, M, step_size, current_timestep, phi, phi_boundary);
        

        /**** visualization code  ****/

        std::cout << "N_dofs " << dofh.NumDofs() << std::endl;


        auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
        for (int global_idx = 0; global_idx < next_timestep.rows(); global_idx++) {
            nodal_data->operator()(dofh.Entity(global_idx)) = next_timestep[global_idx];
        }
        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        vtk_writer.setBinary(true);
        vtk_writer.WritePointData("A", *nodal_data);

        std::cout << "solution norm : " << next_timestep.norm() << std::endl; 

        Eigen::VectorXd backwards_difference  = (next_timestep - current_timestep) / step_size;
        std::cout <<"Backwards difference: " <<  backwards_difference.lpNorm<Eigen::Infinity>() << std::endl; 

        lf::mesh::utils::CodimMeshDataSet<double> relative_permeability{mesh_p, 0, -1};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            int material_tag = cell_tag(*cell); 
            auto material = MaterialFactory::Create(material_tag); 

            Eigen::Vector2d center_of_triangle;
             center_of_triangle << 0.5 , 0.5; 
            auto magnetic_flux = mf_curl(*cell, center_of_triangle);
            Eigen::VectorXd B_field = magnetic_flux[0];
            double reluctivity = material->getReluctivity(B_field.lpNorm<Eigen::Infinity>());
            relative_permeability(*cell) = (1 / reluctivity) / 1.256e-6;
        }

        lf::fe::MeshFunctionGradFE<double, double> mf_grad_solution(fe_space,next_timestep);
        utils::MeshFunctionCurl2DFE mf_curl_solution(mf_grad_solution);

        Eigen::Vector2d init_value = Eigen::Vector2d::Zero();
        lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> H{mesh_p, 0, init_value};
        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
            int material_tag = cell_tag(*cell); 
            auto material = MaterialFactory::Create(material_tag); 

            Eigen::Vector2d center_of_triangle;
             center_of_triangle << 0.5 , 0.5; 
            auto magnetic_flux = mf_curl(*cell, center_of_triangle);
            Eigen::VectorXd B_field = magnetic_flux[0];
            double reluctivity = material->getReluctivity(B_field.lpNorm<Eigen::Infinity>());

            H(*cell) = B_field * reluctivity; 
        }

        lf::fe::MeshFunctionFE<double, double> mf_backwards_difference(fe_space, backwards_difference); 
        lf::mesh::utils::CodimMeshDataSet<double> induced_current{mesh_p, 0, -1};

        for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {

            Eigen::Vector2d center_of_triangle;
            center_of_triangle << 0.5 , 0.5; 
            induced_current(*cell) = - mf_backwards_difference(*cell, center_of_triangle)[0] * cell_conductivity(*cell);
        }

        vtk_writer.WriteCellData("Material", cell_tag);
        vtk_writer.WriteCellData("Backwards-difference", mf_backwards_difference);

        vtk_writer.WriteCellData("B_current", mf_curl_solution);
        vtk_writer.WriteCellData("Conductivity", cell_conductivity);
        vtk_writer.WriteCellData("Prescribed_Current", cell_current);
        vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
        vtk_writer.WriteCellData("H_field", H);

        current_timestep = next_timestep;
        std::cout << std::endl;
    }
    return 0; 
}