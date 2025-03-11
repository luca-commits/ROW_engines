#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"
#include "bdf2.hpp"
#include <fstream>

#include <fstream>

int main (int argc, char *argv[]){

    std::ifstream infile("xinput.txt");
    std::string line;
    double step_size, total_time, max_current, conductivity_ring, exitation_current_parameter;
    std::string mesh_name, exitation_current_type, timestepping_method;
    // when benchmark_creation is 1, the timestep size will be set to 1/10 of the origial size, 
    // and only 1/10 of the visualizations/benchmark files will be saved, so that they can be compared
    // with the original method
    bool benchmark_mode; 
    double benchmark_original_ratio = 2; 
    
    if (!infile.is_open()) {
        std::cerr << "Unable to open the file!" << std::endl;
        return 1;
    }

    // Function to read next non-comment line
    auto readNextValue = [&](auto& var) {
        while (std::getline(infile, line)) {
            if (!line.empty() && line[0] != '#') {
                std::istringstream iss(line);
                if (!(iss >> var)) {
                    std::cerr << "Error reading value: " << line << std::endl;
                    return false;
                }
                return true;
            }
        }
        return false;
    };

    // Read values in correct order
    if (!readNextValue(step_size)) return 1;
    if (!readNextValue(total_time)) return 1;
    if (!readNextValue(mesh_name)) return 1;
    if (!readNextValue(max_current)) return 1;
    if (!readNextValue(conductivity_ring)) return 1;
    if (!readNextValue(exitation_current_type)) return 1;
    if (!readNextValue(exitation_current_parameter)) return 1;
    if (!readNextValue(timestepping_method)) return 1;
    if (!readNextValue(benchmark_mode)) return 1;

    infile.close();

    // Print read values
    std::cout << "Step Size: " << step_size << "\n";
    std::cout << "Total Time: " << total_time << "\n";
    std::cout << "Mesh Name: " << mesh_name << "\n";
    std::cout << "Max Current: " << max_current << "\n";
    std::cout << "Conductivity: " << conductivity_ring << "\n";
    std::cout << "Excitation Type: " << exitation_current_type << "\n";
    std::cout << "Excitation Parameter: " << exitation_current_parameter << "\n";
    std::cout << "Timestepping Method: " << timestepping_method << "\n";
    std::cout << "Benchmark mode" << benchmark_mode << std::endl; 

    if(benchmark_mode){
        step_size /= benchmark_original_ratio; 
    }

    std::cout << "Step Size: " << step_size << "\n";
    std::string benchmark_filename;
    if (benchmark_mode){
        benchmark_filename = std::string("benchmark_file_newton_methods_small_timetep_") + std::string(mesh_name) + std::string(".csv");
    }
    else{
        benchmark_filename = std::string("benchmark_file_newton_methods_") + std::string(mesh_name) + std::string(".csv");
    }

    std::ofstream benchmark_file(benchmark_filename, std::ios::trunc);
    benchmark_file << "# method: " << timestepping_method << ", step size: " << step_size << std::endl; 
    benchmark_file.close();

    bool bdf2 = timestepping_method == "bdf_2";

    std::cout << std::endl << "timestep: " << step_size << std::endl;

    std::filesystem::path here = __FILE__;
    std::string type_of_rotor = mesh_name;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor + "/";
    std::cout << mesh_path;

    auto time_to_current = [max_current, exitation_current_parameter, exitation_current_type](double time){
        if (exitation_current_type == "ramp_up"){
            double rate = 1/exitation_current_parameter;
            return time < exitation_current_parameter ? rate * time * max_current : max_current;
        }
        else{
            const double PI = 3.141592653589793;
            return max_current * std::sin(2 * PI * exitation_current_parameter * time);
        }
    };

    double rel_angle = 0;
    remeshAirgap(mesh_path + "airgap.geo",  mesh_path + "airgap.msh",  M_PI * rel_angle);
    rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
    std::string final_mesh = mesh_path + "motor_" + std::to_string(0) + ".msh";
    unsigned numStableNodes = mergeEverything(mesh_path +  "rotor.msh", mesh_path + "stator.msh", mesh_path + "airgap.msh", final_mesh); 

    std::map<int, double>      tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap
    // std::map<int, double> tag_to_permeability{{1,1.00000037 * MU_0}, {2,0.999994 * MU_0}, {3, 0.999994 * MU_0}, {4,1.00000037 * MU_0}};
    std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}, {4, 0}};
    std::map<int, double> tag_to_conductivity_precoditioner{{1, 400*1e-5}, {2, 400*1e-5}, {3, conductivity_ring}, {4, 400*1e-5}};

    auto [mesh_p_temp, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);


    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());

    std::size_t number_stable_dofs = utils::computeDofsWithoutAirgap(mesh_p_temp, cell_tag, fe_space_temp);


    LF_ASSERT_MSG((number_stable_dofs == numStableNodes), "Something went wrong, the number of stable dofs does not correspond to the number of stableNodes");

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(number_stable_dofs); 
    Eigen::VectorXd previous_timestep = Eigen::VectorXd::Zero(number_stable_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    double angle_step = 2 * M_PI / 60;  //360 / 60 = 6 degrees per timestep

    double time = 0; 
  
    for (unsigned i = 1; time <= total_time; ++i, time += step_size){

        double newton_residual = 1000; 


        double rel_angle = 0; //angle_step * i;
        std::cout << "angle : " << rel_angle << std::endl;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh", rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", rel_angle);
        std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
        unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
 
        std::cout << "timestep: " << i << std::endl;
        std::string vtk_method_name = bdf2 ? "/bdf2" : "/bdf1";
        std::string vtk_filename ;
        if (benchmark_mode){
             vtk_filename = std::string("vtk_files/time_dependent/non-linear/" +   mesh_name  + vtk_method_name + "/eddy_solution_transient_small_timestep") + "_" + std::to_string(i) + std::string(".vtk");
        }else{
            vtk_filename = std::string("vtk_files/time_dependent/non-linear/" +   mesh_name  + vtk_method_name + "/eddy_solution_transient") + "_" + std::to_string(i) + std::string(".vtk");
        }
        double time = i * step_size; 
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 

        auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
        std::cout << "N dofs: " << dofh.NumDofs() << std::endl;

        Eigen::VectorXd current_timestep_extended =  Eigen::VectorXd::Zero(dofh.NumDofs()); 
        Eigen::VectorXd previous_timestep_extended = Eigen::VectorXd::Zero(dofh.NumDofs()); 

        std::cout << "current_timestep.norm() " << current_timestep.norm() << std::endl;  
        current_timestep_extended.head(number_stable_dofs) = current_timestep;       
        lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_timestep_extended);
        utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);
        

        lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_preconditioner{mesh_p, 0, -1};
        // Fill preconditioner conductivity data using the original reader's physical entity numbers
        for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
            // Use the physical entity numbers from the original mesh reading
            unsigned phys_nr = cell_tag(*cell);  // Since we stored this in cell_tag earlier
            cell_conductivity_preconditioner(*cell) = tag_to_conductivity_precoditioner[phys_nr];
        }
            
        Eigen::VectorXd next_newton_step;
        Eigen::VectorXd current_newton_step = current_timestep_extended; //extended          
        double newton_tolerance = 1e-6; 


        for (unsigned j = 1; newton_residual >  newton_tolerance; ++j){ //j < newton_steps

            std::cout << "Newton step " << j << std::endl;
            lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_newton_step);
            utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

            auto [N, rho] = eddycurrent::N_rho_assembler(mesh_p, cell_tag, mf_curl, mf_grad);
            // std::cout << "current newton step - current timestep norm : " << (current_newton_step - current_timestep).norm() << std::endl;

            lf::fe::MeshFunctionGradFE<double, double> mf_grad_current_newton_step(fe_space, current_newton_step); //extended
            utils::MeshFunctionCurl2DFE mf_curl_current_newton_step(mf_grad_current_newton_step);

            auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_current_newton_step);
            auto [A_preconditioner, M_preconditioner, phi_preconditioner] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity_preconditioner, cell_tag, mf_curl_temp);

            if (i == 1 || !bdf2 ) {
                next_newton_step = current_newton_step - newton_step(N, A, M, step_size, current_newton_step,  current_timestep_extended, phi, M_preconditioner); //  current_newton_step
            }
            else next_newton_step = current_newton_step -           bdf2::newton_step(N, A, M, step_size, current_newton_step, previous_timestep_extended, current_timestep_extended, phi, M_preconditioner);

            lf::fe::MeshFunctionGradFE<double, double> mf_grad_next(fe_space, next_newton_step);
            utils::MeshFunctionCurl2DFE mf_curl_next(mf_grad_next);

            auto [A_next, M_next, phi_next] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_next);


            if (i == 1 || !bdf2 ){
                Eigen::SparseMatrix<double> lhs = (step_size * A_next + M);
                Eigen::VectorXd rhs =  M * current_timestep_extended + step_size * phi; //extended
                newton_residual = (lhs * next_newton_step - rhs).norm() / rhs.norm(); 
                std::cout << "newton residual " << newton_residual << std::endl;
 
            }
            else{
                Eigen::SparseMatrix<double> lhs = M + 2./3. * step_size * (A_next);
                Eigen::VectorXd rhs =  + 4./3. * M * current_timestep_extended - 1./3. * M * previous_timestep_extended + 2./3. * step_size * phi;
                newton_residual = (lhs * next_newton_step - rhs).norm() / rhs.norm(); 
                std::cout << "newton residual " << newton_residual << std::endl;
            }

            current_newton_step = next_newton_step; 
            std::cout << std::endl;

        }

        Eigen::VectorXd next_timestep = current_newton_step;

        lf::fe::MeshFunctionGradFE<double, double> mf_grad_solution(fe_space, next_timestep);
        utils::MeshFunctionCurl2DFE mf_curl_solution(mf_grad_solution);


        auto lhs = (step_size * A + M);
        auto rhs =  M * current_timestep_extended + step_size * phi; //extended
        double rel_residual = (lhs * next_timestep - rhs).norm()  / rhs.norm(); 

        if (rhs.norm() > 1e-15) std::cout << "Final newton residual " << rel_residual << std::endl; 


        lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
        vtk_writer.setBinary(true);
        Eigen::VectorXd discrete_solution = next_timestep;


        lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, discrete_solution);
        utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

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
        
        Eigen::VectorXd backwards_difference = (next_timestep  - current_timestep_extended) / step_size;

        std::cout << "next timestep norm : " << next_timestep.norm() << std::endl;
        std::cout << "current timestep norm : " << current_timestep_extended.norm() << std::endl;
        std::cout << "Backwards difference norm : " << backwards_difference.norm() << std::endl;


        //visualize current timestep norm 
        if(i % int(benchmark_original_ratio) == 0 || benchmark_mode == 0){
            std::cout << "writing out visualisation file" << std::endl; 

            auto current_timestep_mesh = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
            for (int global_idx = 0; global_idx < discrete_solution.rows(); global_idx++) {
                current_timestep_mesh->operator()(dofh.Entity(global_idx)) = current_timestep_extended[global_idx];
            }
            vtk_writer.WritePointData("previous_timestep_A", *current_timestep_mesh);


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

            vtk_writer.WriteCellData("induced-current", induced_current);
            std::cout <<"Backwards difference: " <<  backwards_difference.lpNorm<Eigen::Infinity>() << std::endl; 

            //visualize backward difference
            auto backwards_difference_mesh = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
            for (int global_idx = 0; global_idx < discrete_solution.rows(); global_idx++) {
                backwards_difference_mesh->operator()(dofh.Entity(global_idx)) = backwards_difference[global_idx];
            }

            
            benchmark_file.open(benchmark_filename, std::ios::app);
            // benchmark_file << ",";

            // ok now a part that is a bit sketchy but should work. I want to write the power loss to a file, and compare the power
            // losses of two methods. Since the current is a cell based quantity, I will use the cell numbering, and assume it will 
            // stay the same between simulations. Actually this isn't sketchy, because why should the numbering change between
            // simulations (when using same mesh)? 
            for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
                Eigen::Vector2d center_of_triangle;
                center_of_triangle << 0.5 , 0.5;
                benchmark_file << std::pow(mf_backwards_difference(*cell, center_of_triangle)[0], 2) * cell_conductivity(*cell) << "," ; 
            }

            benchmark_file << std::endl;
            benchmark_file.close();

            vtk_writer.WritePointData("Backwards_difference", *backwards_difference_mesh);
            vtk_writer.WriteCellData("B", mf_curl);
            vtk_writer.WriteCellData("Conductivity", cell_conductivity);
            vtk_writer.WriteCellData("Prescribed_Current", cell_current);
            vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
            vtk_writer.WriteCellData("H_field", H);

        }

        std::cout << "current_timestep.size() " << current_timestep.size() << std::endl;
        std::cout << "next_timestep.size() " << next_timestep.size() << std::endl;

        previous_timestep = current_timestep_extended.head(number_stable_dofs);
        current_timestep = next_timestep.head(number_stable_dofs);
        std::cout << std::endl; 

    }

    return 0; 
}