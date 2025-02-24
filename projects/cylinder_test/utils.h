#ifndef UTILS_H_
#define UTILS_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/mesh/mesh.h>
#include <vector>
#include <unordered_map>

namespace utils{

template<class SCALAR_FE, class SCALAR_COEFF> 
class MeshFunctionCurl2DFE{
    public:
        using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));
        MeshFunctionCurl2DFE(lf::fe::MeshFunctionGradFE<SCALAR_FE, SCALAR_COEFF > gradients) : mf_gradients_(gradients){}
        MeshFunctionCurl2DFE(const MeshFunctionCurl2DFE&) = default; 

        std::vector<Eigen::Matrix< Scalar, Eigen::Dynamic, 1>> operator() (const lf::mesh::Entity &e, const Eigen::MatrixXd &local) const {
            std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> result(local.cols());
            const std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> gradients = mf_gradients_(e, local);
            for (Eigen::Index i = 0; i < local.cols(); ++i){
                Eigen::Vector2d temp;
                temp[0] = gradients[i][1];
                temp[1] = -gradients[i][0];
                result[i] = temp;
                // std::cout << "curl norm " << temp.norm() << std::endl; 
            }
            return result; 
        }
    private:
        const lf::fe::MeshFunctionGradFE<SCALAR_FE, SCALAR_COEFF > mf_gradients_;
};

template<class SCALAR_FE, class SCALAR_COEFF> 
class MeshFunctionH{
    public:
        using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));
        MeshFunctionH(const MeshFunctionCurl2DFE<SCALAR_FE, SCALAR_COEFF> curl, lf::mesh::utils::CodimMeshDataSet<SCALAR_COEFF> permeability) : mf_curl_(curl), permeability_(permeability) {}
        std::vector<Eigen::Matrix< Scalar, Eigen::Dynamic, 1>> operator() (const lf::mesh::Entity &e, const Eigen::MatrixXd &local) const {
            auto result = mf_curl_(e, local);
            double permeability = permeability_(e);
            // Divide each matrix in the vector by the scalar
            for (auto& matrix : result) {
                if (permeability/(1.26e-6) > 100){
                    std::cout << "matrix : " << std::endl  << matrix << std::endl;
                    std::cout << "permeability : " << permeability << std::endl; 
                }
                matrix /= permeability;
            }
            return result;
        }
    private:
        const MeshFunctionCurl2DFE<SCALAR_FE, SCALAR_COEFF > mf_curl_;
        lf::mesh::utils::CodimMeshDataSet<double> permeability_;
};

inline size_t computeDofsWithoutAirgap(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    lf::mesh::utils::CodimMeshDataSet<unsigned> cell_tag,
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space) {
    
    // Get the DOF handler
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    
    // Map to store nodes and all their adjacent cell tags
    std::map<const lf::mesh::Entity*, std::set<unsigned>> node_adjacent_tags;
    
    // First pass: collect all adjacent cell tags for each node
    for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
        unsigned tag = cell_tag(*cell);
        auto nodes = cell->SubEntities(2);
        for (const lf::mesh::Entity* node : nodes) {
            node_adjacent_tags[node].insert(tag);
        }
    }
    
    // Set to store the DOFs we want to count
    std::set<lf::base::size_type> stable_dofs;
    
    // Second pass: collect DOFs excluding pure ring nodes
    for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
        unsigned tag = cell_tag(*cell);
        
        // Process nodes for non-ring elements
        if (tag != 4) {
            auto nodes = cell->SubEntities(2);
            for (const lf::mesh::Entity* node : nodes) {
                // Include the node if either:
                // 1. It's a pure non-ring node, or
                // 2. It's on a boundary (has multiple tags)
                if (node_adjacent_tags[node].count(4) == 0 || 
                    node_adjacent_tags[node].size() > 1) {
                    stable_dofs.insert(dofh.GlobalDofIndices(*node).begin(),
                                     dofh.GlobalDofIndices(*node).end());
                }
            }
        }
    }
    
    return stable_dofs.size();
}



// this function computes the area of the sum of all cells where cell_tag == 1 
inline double computeArea(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    const lf::mesh::utils::CodimMeshDataSet<bool>& cell_tag){
        double total_area = 0.0;
        for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
            if (cell_tag(*cell)) {
                const lf::geometry::Geometry* geo_ptr = cell->Geometry();
                total_area += lf::geometry::Volume(*geo_ptr);
            }
        }
        return total_area;
}



class EdgeToCellFinder {
private:
    // Hash function for storing edge endpoints as key
    struct EdgeHasher {
        std::size_t operator()(const std::pair<size_t, size_t>& edge) const {
            return std::hash<size_t>{}(edge.first) ^ 
                   (std::hash<size_t>{}(edge.second) << 1);
        }
    };
    
    // Map to store edge -> cell associations
    std::unordered_map<std::pair<size_t, size_t>, 
                      std::vector<const lf::mesh::Entity*>, 
                      EdgeHasher> edge_to_cell_map;

public:
    explicit EdgeToCellFinder(const lf::mesh::Mesh& mesh) {
        for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
            for (const lf::mesh::Entity* edge : cell->SubEntities(1)) {
                std::span<const lf::mesh::Entity *const> edge_nodes = edge->SubEntities(1);
                size_t v0, v1;
                unsigned int counter = 0;
                for (const lf::mesh::Entity* node : edge_nodes) {
                    if (counter == 0){
                        v0 = mesh.Index(*node);
                    } else {
                        v1 = mesh.Index(*node);
                    }
                    counter++;
                }
                    
                // Store edge endpoints in sorted order for consistent hashing
                auto edge_key = v0 < v1 ? 
                    std::make_pair(v0, v1) : 
                    std::make_pair(v1, v0);
                
                // Add cell to the list of cells containing this edge
                edge_to_cell_map[edge_key].push_back(cell);
            }
        }
    }


    // Find cells adjacent to an edge in O(1) time
    std::vector<const lf::mesh::Entity*> findAdjacentCells(
        const lf::mesh::Entity& edge, 
        const lf::mesh::Mesh& mesh) const {
        
        // Get edge endpoints
        auto edge_nodes = edge.SubEntities(1);
        size_t v0 = mesh.Index(*edge_nodes[0]);
        size_t v1 = mesh.Index(*edge_nodes[1]);
        
        // Create edge key in sorted order
        auto edge_key = v0 < v1 ? 
            std::make_pair(v0, v1) : 
            std::make_pair(v1, v0);
        
        // Look up cells in hash map
        auto it = edge_to_cell_map.find(edge_key);
        if (it != edge_to_cell_map.end()) {
            return it->second;
        }
        return {};
    }
};

// //function that given mesh_p and cell_tag returns two sets of dofs: one corresponding to conducting cells (conductivity > 0) 
// //and one corresponding to non-conducting cells (conductivity == 0)
// //dofs that are both in conducting and non-conducting cells are considered as non-conducting
// inline std::pair<std::set<lf::base::size_type>, std::set<lf::base::size_type>> getConductingDofs(
//     const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
//     const lf::mesh::utils::CodimMeshDataSet<double>& cell_conductivity,
//     const lf::uscalfe::FeSpaceLagrangeO1<double>& fe_space) {
    
//     const lf::assemble::DofHandler& dofh{fe_space.LocGlobMap()};
    
//     std::set<lf::base::size_type> conducting_dofs;
//     std::set<lf::base::size_type> non_conducting_dofs;
    
//     for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
//         double conductivity = cell_conductivity(*cell);
//         auto nodes = cell->SubEntities(2);
//         for (const lf::mesh::Entity* node : nodes) {
//             lf::base::size_type dof = dofh.GlobalDofIndices(*node)[0];
//             //only add the dof to the set if it is not already in the non_conducting_dofs set
//             // if conductivity == 0, check if the dof is already in the conducting_dofs set, if it is, remove it
//             if (conductivity == 0) 
//             {
//                 if (conducting_dofs.contains(dof)) {
//                     conducting_dofs.erase(dof);
//                 }
//                 non_conducting_dofs.insert(dof); 
//             } else {
//                 if (!non_conducting_dofs.contains(dof)) {
//                     conducting_dofs.insert(dof);
//                 }
//             }
//         }
//     }
    
//     return {conducting_dofs, non_conducting_dofs};
// }

} //namespace utils

#endif // UTILS_H