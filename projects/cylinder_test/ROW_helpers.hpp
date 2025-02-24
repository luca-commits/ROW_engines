#include "eddycurrent.h"
#include <eigen3/Dense>

//PRE:
//this function takes a mesh file which contains the mesh file in .gmsh format, the cell current, conductivity, tags, and the magnetic flux density, and three angles for which we want 
// to compute the A, M, and rho matrices
// POST:
// this function returns a tuple of three vectors: the first vector contains three matrices, which are A matrix for three different angles,
// the second vector contains three matrices, which are M matrix for three different angles, the third vector contains three vectors, which are rho vector for three different angles
inline std::tuple<std::vector<Eigen::MatrixXd>, std::vector::Eigen::VectorXd> A_M_rho_assembler(
  std::shared_ptr<const lf::mesh::Mesh> mesh_p,
  lf::mesh::utils::CodimMeshDataSet<double> & cell_conductivity,
  lf::mesh::utils::CodimMeshDataSet<unsigned int> & cell_tags,
  utils::MeshFunctionCurl2DFE<double, double>  & cell_B,
  std::vector<double> angles, 
  std::string mesh_path,
  std::map<int, double> tag_to_conductivity
  ){
    
    std::vector<Eigen::SparseMatrix<double>> A_matrices; 
    std::vector<Eigen::VectorXd> rho_vectors; 

    std::map<int, double> tag_to_current{}; //1 -> air, 2-> cylinder, 3 -> ring, 4 -> airgap

    auto [mesh_p_temp, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);

    for (auto angle : angles){
      remeshAirgap(mesh_path + "airgap.geo", mesh_path + "airgap.msh", rel_angle);
      rotateAllNodes_alt(mesh_path + "rotor.msh", mesh_path + "rotated_rotor.msh", rel_angle);
      std::string final_mesh = mesh_path + "motor_" + std::to_string(i) + ".msh";
      unsigned numStableNodes = mergeEverything(mesh_path +  "stator.msh", mesh_path + "rotated_rotor.msh", mesh_path + "airgap.msh",final_mesh);
      std::cout<< "Final mesh: " << final_mesh << std::endl;

      std::cout << "timestep: " << i << std::endl;
      std::string vtk_filename = std::string("vtk_files/time_dependent/rotating/" + mesh_name  + "/eddy_solution_transient_") + "_" + std::to_string(i) + std::string(".vtk");
      double time = i * step_size; 
      std::cout << "current " << time_to_current(time) << std::endl;
      tag_to_current = {{1,0},  {2, time_to_current(time)}, {3, 0}, {4, 0}}; 

      auto [mesh_p, cell_current, cell_conductivity, cell_tag] = eddycurrent::readMeshWithTags(final_mesh, tag_to_current, tag_to_conductivity);
      auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

      const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
      std::cout << "N dofs: " << dofh.NumDofs() << std::endl;

      Eigen::VectorXd current_timestep_extended = Eigen::VectorXd::Zero(dofh.NumDofs()); 

      std::cout << "current_timestep.norm() " << current_timestep.norm() << std::endl;        
      current_timestep_extended.head(number_stable_dofs) = current_timestep; 
      std::cout << "current_timestep_extended.norm() " << current_timestep_extended.norm() << std::endl;
      lf::fe::MeshFunctionGradFE<double, double> mf_grad_temp(fe_space, current_timestep_extended);
      utils::MeshFunctionCurl2DFE mf_curl_temp(mf_grad_temp);

      A_matrices.push_back(A); 

      std::map<int, double> tag_to_conductivity{{1, 0}, {2, 0}, {3, conductivity_ring}, {4, 0}};

    }
  }