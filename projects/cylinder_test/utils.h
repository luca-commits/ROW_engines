#ifndef UTILS_H_
#define UTILS_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

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
                matrix /= permeability;
            }
            return result;
        }
    private:
        const MeshFunctionCurl2DFE<SCALAR_FE, SCALAR_COEFF > mf_curl_;
        lf::mesh::utils::CodimMeshDataSet<double> permeability_;
};

size_t computeDofsWithoutAirgap(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
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

} //namespace utils

#endif // UTILS_H