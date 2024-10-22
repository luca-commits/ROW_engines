/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements from a Gmsh generated mesh
 * @author Simon Meierhans
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  std::ofstream out("nodes.txt");

  // find path to mesh
  std::filesystem::path here = __FILE__;
  auto mesh_path = here.remove_filename() / "meshes/cylinder_spherical.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());
  auto mesh = reader.mesh();
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  lf::io::VtkWriter vtk_writer(mesh, "vtk_files/smiley.vtk");
  vtk_writer.setBinary(1);
  auto mds = lf::mesh::utils::make_CodimMeshDataSet<int>(mesh, 2, 0);
  // mark the eyes
  for (const auto* e : mesh->Entities(2)) {
    mds->operator()(*e) = 1;
  }
  vtk_writer.WritePointData("smile", *mds);
  out.close();

  lf::io::writeTikZ(*mesh, "smiley.tikz", [](const lf::mesh::Entity & e){return true;});

}
