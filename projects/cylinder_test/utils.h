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
        std::vector<Eigen::Matrix< Scalar, Eigen::Dynamic, 1>> operator() (const lf::mesh::Entity &e, const Eigen::MatrixXd &local) const {
            std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> result(local.cols());
            const std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> gradients = mf_gradients_(e, local);
            for (Eigen::Index i = 0; i < local.cols(); ++i){
                Eigen::Vector2d temp;
                //temp << gradients[i][1] , -gradients[i][0];
                temp[0] = gradients[i][1];
                temp[1] = -gradients[i][0];
                result[i] = temp;
                // if (i == 0) {std::cout << "gradients shape" << gradients[i].size() << std::endl;}

                // result[i][0] = gradients[i][1];
                // result[i][1] = -gradients[i][0];
            }
            return result; 
        }
    private:
        lf::fe::MeshFunctionGradFE<SCALAR_FE, SCALAR_COEFF > mf_gradients_;
};
} //namespace utils

#endif // UTILS_H