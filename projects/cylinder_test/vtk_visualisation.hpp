
lf::io::VtkWriter vtk_writer(mesh_p, vtk_filename);
vtk_writer.setBinary(true);

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

vtk_writer.WriteCellData("B", mf_curl);
vtk_writer.WriteCellData("Conductivity", cell_conductivity);
vtk_writer.WriteCellData("Prescribed_Current", cell_current);
vtk_writer.WriteCellData("Relative_Permeability", relative_permeability);
vtk_writer.WriteCellData("H_field", H);