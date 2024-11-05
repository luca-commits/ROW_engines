#include <Eigen/Dense>
#include <gmsh.h> 
#include <iostream> 


int test_coincidence(std::string mesh_file){

    gmsh::initialize();
        
    // Open the mesh file
    gmsh::merge(mesh_file);
    // Get all nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;
    
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);

    std::cout << "Point 73" << std::endl << "x " << coords[73 * 3] << "y " << coords[73 * 3 + 1];

    bool upsie = 0;
    
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        for (std::size_t j = 0; j < nodeTags.size(); ++j){
            double x_first = coords[i * 3];
            double y_first = coords[i * 3 + 1];
            double x_second = coords[j * 3];
            double y_second = coords[j * 3 + 1];
            
            Eigen::Vector2d first;
            first << x_first, y_first;
            Eigen::Vector2d second; 
            second << x_second, y_second;

            if ((first -second).norm() < 1e-7){
                if (i != j){
                    std::cout << "Node difference " << (first - second).norm() << std::endl;
                    std::cout << "first coordinates " << std::endl << first << std::endl;
                    std::cout << "second coordinates " << std::endl << second << std::endl; 
                    std::cout << "i " << i << ", j " << j << std::endl; 
                    upsie = 1; 
                } 
            }
        }
    }

    if (upsie) return -1; 

    gmsh::finalize(); 
    return 0; 
}

int main(){
    if(test_coincidence("motor_0.msh") != 0){
        std::cout << "Something went wrong" << std::endl;
    }
    else 
        std::cout << "All good so far" << std::endl; 
    return 0; 
}