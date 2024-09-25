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

int main() {
  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  // find path to mesh
  std::filesystem::path here = __FILE__;
  auto mesh_path = here.remove_filename() / "meshes/t30.msh";

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());

  unsigned nr_physical_entities = reader.PhysicalEntities(0).size();

  std::cout << nr_physical_entities << std::endl;

  // get pointer to mesh
  auto mesh = reader.mesh();

  lf::io::VtkWriter vtk_writer(mesh, "t30_mesh_lehrfem.vtk");
  auto cell_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh, 1);
  vtk_writer.WriteCellData("cellData", *cell_data);

}
