#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"
#include "rosenbrock_wanner.hpp"

#include <fstream>

int main (int argc, char *argv[]){


    std::ifstream infile("xinput.txt");
    std::string line;
    double step_size, total_time, max_current, conductivity_ring, exitation_current_parameter, adaptive_rel_tol = 1e-1, adaptive_abs_tol = 0.02, strict_ratio = 5.;
    std::string mesh_name, exitation_current_type, timestepping_method, geometry_type;
    bool adaptive; 
    
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
    if (!readNextValue(geometry_type)) return 1; 
    if (!readNextValue(timestepping_method)) return 1;
    if (!readNextValue(adaptive)) return 1;

    infile.close();

    // Print read values
    std::cout << "Step Size: " << step_size << "\n";
    std::cout << "Total Time: " << total_time << "\n";
    std::cout << "Mesh Name: " << mesh_name << "\n";
    std::cout << "Max Current: " << max_current << "\n";
    std::cout << "Conductivity: " << conductivity_ring << "\n";
    std::cout << "Excitation Type: " << exitation_current_type << "\n";
    std::cout << "Excitation Parameter: " << exitation_current_parameter << "\n";
    std::cout << "Timestepping Method: " << "rosenbrock wanner" << "\n";
    std::cout << "Adaptive:  " << adaptive << "\n"; 

    std::cout << "Geometry type : " << geometry_type << std::endl; 
    double original_step_size = step_size;


    std::string benchmark_filename_current = std::string("benchmark_file_rosenbrock-wanner_") + std::string(mesh_name) + std::string("_current.csv");
    std::string benchmark_filename_B_field = std::string("benchmark_file_rosenbrock-wanner_") + std::string(mesh_name) + std::string("_magnetic_field.csv");
    std::string benchmark_filename_B_power = std::string("benchmark_file_rosenbrock-wanner_") + std::string(mesh_name) + std::string("_magnetic_power.csv");

    std::ofstream benchmark_file_current(benchmark_filename_current, std::ios::trunc);
    benchmark_file_current << "# method: " << "rosenbrock wanner" << ", step size: " << step_size << std::endl; 
    benchmark_file_current.close();
    
    std::ofstream benchmark_file_B_field(benchmark_filename_B_field, std::ios::trunc);
    benchmark_file_B_field << "# method: " << "rosenbrock wanner" << ", step size: " << step_size << std::endl; 
    benchmark_file_B_field.close();

    std::ofstream benchmark_file_B_power(benchmark_filename_B_power, std::ios::trunc);
    benchmark_file_B_power << "# method: " << "rosenbrock wanner" << ", step size: " << step_size << std::endl; 
    benchmark_file_B_power.close();


    bool rotating_geometry = true;
    if (geometry_type == "transformer") rotating_geometry = false; 
    else rotating_geometry = true; 

    std::cout << "rotating? : " << rotating_geometry << std::endl; 

    std::filesystem::path here = __FILE__;
    std::string type_of_rotor = mesh_name;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor + "/";
    std::cout << mesh_path << std::endl;

    auto time_to_current = [max_current, exitation_current_parameter, exitation_current_type](double current_time){
        if (exitation_current_type == "ramp_up"){
            double rate = 1/exitation_current_parameter;
            return current_time < exitation_current_parameter ? rate * current_time * max_current : max_current;
        }
        else{
            const double PI = 3.141592653589793;
            return  max_current * std::sin(2 * PI * exitation_current_parameter * current_time);
        }
    };

    auto time_to_current_derivative = [max_current, exitation_current_parameter, exitation_current_type](double current_time){
        if (exitation_current_type == "ramp_up"){
            double rate = 1/exitation_current_parameter;
            return current_time < exitation_current_parameter ? rate  * max_current : 0;
        }
        else{
            const double PI = 3.141592653589793;
            return max_current  * std::cos(2 * PI * exitation_current_parameter * current_time);
        }
    };

    double rel_angle = 0;
    std::string final_mesh;
    unsigned numStableNodes;
    if(rotating_geometry){
        remeshAirgap(mesh_path + "airgap.geo",  mesh_path + "airgap.msh",  M_PI * rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
        final_mesh = mesh_path + "motor_" + std::to_string(0) + ".msh";
        numStableNodes = mergeEverything(mesh_path +  "rotor.msh", mesh_path + "stator.msh", mesh_path + "airgap.msh", final_mesh); 
    }
    else{
        final_mesh = std::string(here.remove_filename()) + std::string("meshes/transformator/") + mesh_name + std::string(".msh");
    }


    std::map<int, double>      tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap
    std::map<int, double>      tag_to_current_derivative{}; 
    std::map<int, double> tag_to_conductivity;
    std::map<int, double> tag_to_conductivity_precoditioner;
    if (geometry_type == "transformer"){
        tag_to_conductivity = {{6, conductivity_ring}, {2, 0}, {3, 0}, {4, 0}, {5, 0}};
        tag_to_conductivity_precoditioner = {{6, conductivity_ring}, {2, 400*1e-5}, {3, 400*1e-5}, {4, 400*1e-5}, {5, 400*1e-5}};
    }
    else{
        tag_to_conductivity = {{1, 0}, {2, 0}, {3, conductivity_ring}, {4, 0}};
        tag_to_conductivity_precoditioner = {{1, 400*1e-5}, {2, 400*1e-5}, {3, conductivity_ring}, {4, 400*1e-5}};
    }

    auto [mesh_p_temp, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
    auto fe_space_temp = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_temp);

    lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_preconditioner{mesh_p_temp, 0, -1};
    // Fill preconditioner conductivity data using the original reader's physical entity numbers
    for (const lf::mesh::Entity *cell : mesh_p_temp->Entities(0)) {
        // Use the physical entity numbers from the original mesh reading
        unsigned phys_nr = cell_tag(*cell);  // Since we stored this in cell_tag earlier
        cell_conductivity_preconditioner(*cell) = tag_to_conductivity_precoditioner[phys_nr];
    }
            

    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh_temp{fe_space_temp->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh_temp.NumDofs());


    std::size_t number_stable_dofs;
    if (rotating_geometry){
        number_stable_dofs = utils::computeDofsWithoutAirgap(mesh_p_temp, cell_tag, fe_space_temp);
        LF_ASSERT_MSG((number_stable_dofs == numStableNodes), "Something went wrong, the number of stable dofs does not correspond to the number of stableNodes");
    }
    else {
        number_stable_dofs = N_dofs; 
    }

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(N_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    double angle_step = 0;  //360 / 60 = 6 degrees per timestep
    double current_time = 0; 
    unsigned i; 
    for (; current_time <= total_time;){

        const double PI = 3.141592653589793;
        double rel_angle = 0 ; // angle_step * i;
         
        // std::cout << "angle : " << rel_angle << std::endl;
               std::string final_mesh;
        unsigned numStableNodes;
        if(rotating_geometry){
            remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh", rel_angle);
            rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", rel_angle);
            final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
            numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
        }
        else{
            final_mesh = std::string(here.remove_filename()) + std::string("meshes/transformator/") + mesh_name + std::string(".msh");
        }

        if (geometry_type == "transformer"){
            tag_to_current = {{6,0},  {2, time_to_current(current_time)}, {3, 0}, {4, 0}, {5, -time_to_current(current_time)}};
            tag_to_current_derivative = {{6,0},  {2, time_to_current_derivative(current_time)}, {3, 0}, {4, 0}, {5, -time_to_current_derivative(current_time)}};
        }
        else{
            tag_to_current = {{1,0},  {2, time_to_current(current_time)}, {3, 0}, {4, 0}}; 
            tag_to_current_derivative = {{1,0},  {2, time_to_current_derivative(current_time)}, {3, 0}, {4, 0}};
        }
        

        auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
      
        lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_timestep);
        utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);

        lf::mesh::utils::CodimMeshDataSet<double> cell_conductivity_preconditioner{mesh_p, 0, -1};
        // Fill preconditioner conductivity data using the original reader's physical entity numbers
        for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
            // Use the physical entity numbers from the original mesh reading
            unsigned phys_nr = cell_tag(*cell);  // Since we stored this in cell_tag earlier
            cell_conductivity_preconditioner(*cell) = tag_to_conductivity_precoditioner[phys_nr];
        }

        
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);
        auto [A_preconditioner, M_preconditioner, phi_preconditioner] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity_preconditioner, cell_tag, mf_curl_temp);

        // std::cout << "A norm : " << A.norm() << std::endl ; 


        lf::mesh::utils::CodimMeshDataSet<double> cell_current_derivative = eddycurrent::getCellCurrent(mesh_p, tag_to_current_derivative, cell_tag);
        Eigen::VectorXd time_derivative = eddycurrent::phi_assembler(mesh_p, cell_current_derivative);        
        Eigen::SparseMatrix<double> N = eddycurrent::N_assembler(mesh_p, cell_tag, mf_curl_temp, mf_grad_temp);

        Eigen::SparseMatrix<double> jacobian = -(A + N);

        lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tags = cell_tag; 
        std::shared_ptr<const lf::mesh::Mesh> mesh_ps = mesh_p; 

        auto  function_evaluator = [&max_current, &exitation_current_parameter, &exitation_current_type, 
                                    &cell_tags, &mesh_ps, &step_size, &geometry_type, &tag_to_current]( double time, Eigen::VectorXd current_timestep,int i) //i is for debugging purposes only
                                    -> Eigen::VectorXd {


            auto time_to_current = [max_current, exitation_current_parameter, exitation_current_type](double time){
                if (exitation_current_type == "ramp_up"){
                    double rate = 1/exitation_current_parameter;
                    return time < exitation_current_parameter ? rate * time * max_current : max_current;
                }
                else{
                    const double PI = 3.141592653589793;
                    return max_current * std::sin(2 * PI * exitation_current_parameter * time);
                    std::cout << "sine " << std::sin(2 * PI * exitation_current_parameter * time) << std::endl;
                }
            };

            if (geometry_type == "transformer"){
                tag_to_current = {{6,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}, {5, -time_to_current(time)}};
            }
            else{
                tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
            }
            lf::mesh::utils::CodimMeshDataSet<double> cell_current = eddycurrent::getCellCurrent(mesh_ps, tag_to_current, cell_tags);

            auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_ps);

            lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_timestep);
            utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

            Eigen::SparseMatrix<double> A = eddycurrent::A_assembler(mesh_ps, cell_tags, mf_curl); 
            Eigen::VectorXd rho = A * current_timestep; 
            Eigen::VectorXd load = eddycurrent::phi_assembler(mesh_ps, cell_current);

            if (bool debug = 1){
                std::string vtk_filename = std::string("vtk_files/time_dependent/debug_rho_load_") + std::to_string(i) + std::string(".vtk");
                lf::io::VtkWriter vtk_writer(mesh_ps, vtk_filename);

                vtk_writer.setBinary(true);

                auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_ps);
                const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

                auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_ps, 2);
                for (int global_idx = 0; global_idx < rho.rows(); global_idx++) {
                    nodal_data->operator()(dofh.Entity(global_idx)) = step_size * load[global_idx];
                }
                vtk_writer.WritePointData("load", *nodal_data);

                nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_ps, 2);
                for (int global_idx = 0; global_idx < rho.rows(); global_idx++) {
                    nodal_data->operator()(dofh.Entity(global_idx)) = step_size *  rho[global_idx];
                }
                vtk_writer.WritePointData("rho", *nodal_data);
            }

            return load - rho; 

        };

        auto [next_timestep_high, next_timestep_low] = row_step(step_size, current_time, current_timestep, jacobian, M, M_preconditioner , time_derivative, function_evaluator, mesh_p);
  
        double high_low_diff = (next_timestep_high - next_timestep_low).norm(); 
        bool accept = 0; 
        bool strongly_accept = 0; 
        if (high_low_diff < std::max(adaptive_rel_tol * current_timestep.norm(), adaptive_abs_tol)) accept = true; 
        if (high_low_diff < std::max((adaptive_rel_tol / strict_ratio) * current_timestep.norm(), adaptive_abs_tol / strict_ratio)) strongly_accept = true; 
        if (accept){
            current_time += step_size;
            std::cout << "step accepted" << std::endl; 
        }
        std::cout << "current_time : " << current_time << std::endl;
        std::cout << "std::fmod(current_time, original_step_size)  " << std::fmod(current_time, original_step_size)   << std::endl; 

        //I only want to write out the visualization when the time corresponds to a regular interval such that I can compare results
        //with newton-iteration based methods
        if ((std::fmod(current_time, original_step_size) < 1e-12 && accept == true) || !adaptive){
            unsigned i = int(current_time / original_step_size);
            std::string vtk_filename = std::string("vtk_files/time_dependent/row_static/" + mesh_name) + "_" + std::to_string(i) + std::string(".vtk");
            std::cout << "Writing visulasation file to"<< vtk_filename << std::endl; 
            lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);

            vtk_writer.setBinary(true);
            Eigen::VectorXd discrete_solution = next_timestep_high;

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
            
            Eigen::VectorXd backwards_difference = (next_timestep_high - current_timestep) / step_size;

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


            benchmark_file_current.open(benchmark_filename_current, std::ios::app);
            benchmark_file_B_field.open(benchmark_filename_B_field, std::ios::app);
            benchmark_file_B_power.open(benchmark_filename_B_power, std::ios::app);

            // ok now a part that is a bit sketchy but should work. I want to write the power loss to a file, and compare the power
            // losses of two methods. Since the current is a cell based quantity, I will use the cell numbering, and assume it will 
            // stay the same between simulations. Actually this isn't sketchy, because why should the numbering change between
            // simulations (when using same mesh)? 
            for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {

                int material_tag = cell_tag(*cell); 
                auto material = MaterialFactory::Create(material_tag); 
                Eigen::Vector2d center_of_triangle;
                center_of_triangle << 0.5 , 0.5; 
                auto magnetic_flux = mf_curl(*cell, center_of_triangle);
                Eigen::VectorXd B_field = magnetic_flux[0];

                auto magnetic_flux_previous = mf_curl_previous(*cell, center_of_triangle); 
                Eigen::VectorXd B_field_previous = magnetic_flux_previous[0];

                Eigen::VectorXd B_dot = (B_field - B_field_previous) / step_size; 


                double reluctivity = material->getReluctivity(B_field.lpNorm<Eigen::Infinity>());

                benchmark_file_current << std::pow(mf_backwards_difference(*cell, center_of_triangle)[0], 2) * cell_conductivity(*cell) << "," ; 
                benchmark_file_B_field << 0.5 * B_field.dot(B_field )* reluctivity << "," ; 
                benchmark_file_B_power << reluctivity * B_field.dot(B_dot) + std::pow(mf_backwards_difference(*cell, center_of_triangle)[0], 2) * cell_conductivity(*cell) << "," ; 
            }

            benchmark_file_current << std::endl;
            benchmark_file_current.close();

            benchmark_file_B_field << std::endl;
            benchmark_file_B_field.close();

            benchmark_file_B_power << std::endl;
            benchmark_file_B_power.close();


            vtk_writer.WriteCellData("induced-current", induced_current);
            std::cout <<"Backwards difference: " <<  backwards_difference.lpNorm<Eigen::Infinity>() << std::endl; 

            vtk_writer.WriteCellData("B", mf_curl);
            vtk_writer.WriteCellData("Conductivity", cell_conductivity);
            vtk_writer.WriteCellData("Prescribed_Current", cell_current);
            vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
            vtk_writer.WriteCellData("H_field", H);
        }

        
        std::cout << "high_low_diff: " << high_low_diff << std::endl; 
        std::cout << "relative tolerance: " << (adaptive_rel_tol / strict_ratio) * current_timestep.norm() << std::endl; 
        if (adaptive){
            if (strongly_accept && step_size < original_step_size && std::fmod(current_time, 2 * step_size) < 1e-12){
                step_size *= 2;
                current_timestep = next_timestep_high; 
                std::cout << "step strongly accepted " << std::endl; 
            }
            else if (accept){
                current_timestep = next_timestep_high;  
            }
            else{
                std::cout << "step rejected " << std::endl;
                step_size /= 2; 
                std::cout << "step size reduced to " << step_size << std::endl; 
            }
        }
        else{
            current_timestep = next_timestep_high; 
            current_time += step_size; 
        }
    }
    return 0; 
}