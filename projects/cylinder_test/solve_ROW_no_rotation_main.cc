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
    double step_size, total_time, max_current, conductivity_ring, exitation_current_parameter;
    std::string mesh_name, exitation_current_type, timestepping_method;
    
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


     std::string benchmark_filename = std::string("benchmark_file_rosenbrock-wanner_") + std::string(mesh_name) + std::string(".csv");

    std::ofstream benchmark_file(benchmark_filename, std::ios::trunc);
    benchmark_file << "# method: " << "rosenbrock wanner" << ", step size: " << step_size << std::endl; 
    benchmark_file.close();

    std::filesystem::path here = __FILE__;
    std::string type_of_rotor = mesh_name;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/rotating/") + type_of_rotor + "/";
    std::cout << mesh_path << std::endl;

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

    auto time_to_current_derivative = [max_current, exitation_current_parameter, exitation_current_type](double time){
        if (exitation_current_type == "ramp_up"){
            double rate = 1/exitation_current_parameter;
            return time < exitation_current_parameter ? rate  * max_current : 0;
        }
        else{
            const double PI = 3.141592653589793;
            return max_current * std::cos(2 * PI * exitation_current_parameter * time);
        }
    };



    double rel_angle = 0;
    remeshAirgap(mesh_path + "airgap.geo",  mesh_path + "airgap.msh",  M_PI * rel_angle);
    rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", M_PI * rel_angle);
    std::string final_mesh = mesh_path + "motor_" + std::to_string(0) + ".msh";
    std::cout << "Mesh Path :" << final_mesh << std::endl;
    unsigned numStableNodes = mergeEverything(mesh_path +  "rotor.msh", mesh_path + "stator.msh", mesh_path + "airgap.msh", final_mesh); 

    std::map<int, double>      tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap
    std::map<int, double>      tag_to_current_derivative{}; 
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

    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(N_dofs); 
    std::cout << "Conductivity ring: " << conductivity_ring << std::endl;

    double angle_step = 0;  //360 / 60 = 6 degrees per timestep
    double time = 0; 

    for (unsigned i = 1; time <= total_time; ++i, time += step_size){
        double rel_angle = 0 ; // angle_step * i;
        // std::cout << "angle : " << rel_angle << std::endl;
        remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh", rel_angle);
        rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", rel_angle);
        std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
        unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
        std::cout << "timestep: " << i << std::endl;
        std::string vtk_filename = std::string("vtk_files/time_dependent/row_static/" + mesh_name) + "_" + std::to_string(i) + std::string(".vtk");
        double time = i * step_size; 

        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
        tag_to_current_derivative = {{1,0},  {2, time_to_current_derivative(time)}, {3, 0}, {4, 0}};

        auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
        auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

        const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
      
        lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_timestep);
        utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);
        
        auto [A, M, phi] = eddycurrent::A_M_phi_assembler(mesh_p, cell_current, cell_conductivity, cell_tag, mf_curl_temp);

        lf::mesh::utils::CodimMeshDataSet<double> cell_current_derivative = eddycurrent::getCellCurrent(mesh_p, tag_to_current_derivative, cell_tag);
        Eigen::VectorXd time_derivative = eddycurrent::phi_assembler(mesh_p, cell_current_derivative);        
        Eigen::SparseMatrix<double> N = eddycurrent::N_assembler(mesh_p, cell_tag, mf_curl_temp, mf_grad_temp);

        Eigen::SparseMatrix<double> jacobian = -(A + N);

        lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tags = cell_tag; 
        std::shared_ptr<const lf::mesh::Mesh> mesh_ps = mesh_p; 

        auto  function_evaluator = [&max_current, &exitation_current_parameter, &exitation_current_type, &cell_tags, &mesh_ps, &step_size]( double time, Eigen::VectorXd current_timestep)
                                    -> Eigen::VectorXd {


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

            std::map<int, double> tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 
            lf::mesh::utils::CodimMeshDataSet<double> cell_current = eddycurrent::getCellCurrent(mesh_ps, tag_to_current, cell_tags);

            auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_ps);

            lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, current_timestep);
            utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

            Eigen::VectorXd rho =  eddycurrent::rho_assembler(mesh_ps, cell_tags, mf_curl);
            Eigen::VectorXd load = eddycurrent::phi_assembler(mesh_ps, cell_current);

            return load - rho; 

        };

        Eigen::VectorXd next_timestep = row_step(step_size, time, current_timestep, jacobian, M , time_derivative, function_evaluator);

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
        
        Eigen::VectorXd backwards_difference = (next_timestep - current_timestep) / step_size;

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

        vtk_writer.WriteCellData("induced-current", induced_current);
        std::cout <<"Backwards difference: " <<  backwards_difference.lpNorm<Eigen::Infinity>() << std::endl; 

        vtk_writer.WriteCellData("B", mf_curl);
        vtk_writer.WriteCellData("Conductivity", cell_conductivity);
        vtk_writer.WriteCellData("Prescribed_Current", cell_current);
        vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
        vtk_writer.WriteCellData("H_field", H);

        current_timestep = next_timestep; 
    }
    return 0; 
}