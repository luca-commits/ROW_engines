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


} //namespace utils

#endif // UTILS_H