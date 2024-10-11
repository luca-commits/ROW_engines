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
        infile.close();
    }
    else{
        std::cerr << "Unable to open the file!" << std::endl;
    }

    double step_size = total_time / double(timesteps); 
    std::cout << "timestep: " << step_size << std::endl;

    std::filesystem::path here = __FILE__;
    auto mesh_path = std::string(here.remove_filename()) + std::string("meshes/") + mesh_name;
    std::cout << mesh_path;

    auto time_to_current = [max_current](double time){
        return time < 0.001 ? 1000 * time * max_current : max_current;
    };

    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    // std::cout << "Reading mesh from " << mesh_path << std::endl;
    lf::io::GmshReader reader(std::move(mesh_factory), mesh_path);
    std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};

    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    // Dimension of finite element space = number of nodes of the mesh
    const lf::base::size_type N_dofs(dofh.NumDofs());
    // std::cout << "Solving BVP with " << N_dofs << " FE dofs" << std::endl;


    Eigen::VectorXd current_timestep = Eigen::VectorXd::Zero(N_dofs); 

    std::map<int, double>      tag_to_current{};
    std::map<int, double> tag_to_conductivity{{1,0}, {2,0}, {3, 1e7 }}; //1 -> air, 2-> cylinder, 3 -> ring
    std::map<int, double> tag_to_permeability{{1,1},  {2,1}, {3, 20}};

    for (unsigned i = 0; i < timesteps; ++i){
        std::cout << "timestep: " << i << std::endl;
        double time = i * step_size; 
        std::cout << "current " << time_to_current(time) << std::endl;
        tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}}; //1 -> air, 2-> cylinder, 3 -> ring
        auto [mesh_p, cell_current, cell_permeability, cell_conductivity] = eddycurrent_time_dependent::readMeshWithTags(mesh_path, tag_to_current, tag_to_permeability, tag_to_conductivity);
        auto [A, M, phi] = eddycurrent_time_dependent::A_M_phi_assembler(mesh_p, cell_current, cell_permeability, cell_conductivity);
        // std::cout << "M " << M << std::endl;
        Eigen::VectorXd next_timestep = implicit_euler_step(A, M, step_size, current_timestep, phi);
        current_timestep = next_timestep;
    }

    Eigen::VectorXd discrete_solution = current_timestep; 

    std::string vtk_filename = std::string("vtk_files/time_dependent/eddy_solution_") + argv[1] + std::string(".vtk");

    lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);

    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < discrete_solution.rows();
        global_idx++) {
        nodal_data->operator()(dofh.Entity(global_idx)) =
            discrete_solution[global_idx];
    };


    lf::mesh::utils::CodimMeshDataSet<double> non_zero_current{mesh_p, 0, -1};
    for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
        LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                    " edge must be a triangle!");
        non_zero_current(*cell) = tag_to_current[reader.PhysicalEntityNr(*cell)[0]] != 0;
    }

    lf::mesh::utils::CodimMeshDataSet<double> conductivity{mesh_p, 0, -1};
    for (const lf::mesh::Entity *cell : mesh_p -> Entities(0)) {
        LF_ASSERT_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                    " edge must be a triangle!");
        conductivity(*cell) = tag_to_conductivity[reader.PhysicalEntityNr(*cell)[0]];
    }

    vtk_writer.WritePointData("A*", *nodal_data);
    vtk_writer.WriteCellData("non_zero_current", non_zero_current);
    // vtk_writer.WriteCellData("conductivity", conductivity);

    lf::fe::MeshFunctionGradFE<double, double> mf_grad(fe_space, discrete_solution);
    std::cout << "computing B from A*" << std::endl;
    utils::MeshFunctionCurl2DFE mf_curl(mf_grad);

    vtk_writer.WriteCellData("B", mf_curl);
    return 0; 
}