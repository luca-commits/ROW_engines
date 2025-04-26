#include "implicit_euler.h"
#include "eddycurrent.h"
#include "utils.h"
#include "remesh_airgap.hpp"
#include "rotate_mesh.hpp"

#include <fstream>

int main (int argc, char *argv[]){

    double conductivity = 1;

   Eigen::MatrixXd V (2, 3); 
   V << 0, 0, 1,
        0, 1, 0;
    // hard-coded quadrature rule: 






 // Define the type alias at the beginning
using Vec2d = Eigen::Vector2d;  // or Eigen::Matrix<double, 2, 1>

// Higher-order quadrature rule (12 points)
Eigen::MatrixXd quadrature_points(12, 3);
quadrature_points << 
    0.063089014491502228340331602870, 0.063089014491502228340331602870, 0.873821971016995543319336794260,
    0.063089014491502228340331602870, 0.873821971016995543319336794260, 0.063089014491502228340331602870,
    0.873821971016995543319336794260, 0.063089014491502228340331602870, 0.063089014491502228340331602870,
    0.249286745170910421291638553107, 0.249286745170910421291638553107, 0.501426509658179157416722893786,
    0.249286745170910421291638553107, 0.501426509658179157416722893786, 0.249286745170910421291638553107,
    0.501426509658179157416722893786, 0.249286745170910421291638553107, 0.249286745170910421291638553107,
    0.310352451033784405473666148844, 0.053145049844816947353249671092, 0.636502499121398647173084180064,
    0.053145049844816947353249671092, 0.310352451033784405473666148844, 0.636502499121398647173084180064,
    0.636502499121398647173084180064, 0.310352451033784405473666148844, 0.053145049844816947353249671092,
    0.636502499121398647173084180064, 0.053145049844816947353249671092, 0.310352451033784405473666148844,
    0.310352451033784405473666148844, 0.636502499121398647173084180064, 0.053145049844816947353249671092,
    0.053145049844816947353249671092, 0.636502499121398647173084180064, 0.310352451033784405473666148844;

Eigen::VectorXd quadrature_weights(12);
quadrature_weights << 
    0.050844906370206816921300031891,
    0.050844906370206816921300031891,
    0.050844906370206816921300031891,
    0.116786275726379071434620099573,
    0.116786275726379071434620099573,
    0.116786275726379071434620099573,
    0.082851075618373575193884366397,
    0.082851075618373575193884366397,
    0.082851075618373575193884366397,
    0.082851075618373575193884366397,
    0.082851075618373575193884366397,
    0.082851075618373575193884366397;

// Pre-compute area
const double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) - 
                                 (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));

// Pre-compute the gradients for each basis function
std::vector<Eigen::Vector2d> gradients(3);
gradients[0] = Eigen::Vector2d(V(1,1) - V(1,2), V(0,2) - V(0,1));
gradients[1] = Eigen::Vector2d(V(1,2) - V(1,0), V(0,0) - V(0,2));
gradients[2] = Eigen::Vector2d(V(1,0) - V(1,1), V(0,1) - V(0,0));

const double factor = 1.0 / (2.0 * area);
for(auto& grad : gradients) {
    grad *= factor;
}

// Modified lambda functions using pre-computed gradients
auto lambda1 = [&V, &gradients](const Eigen::Vector2d& x) -> double {
    return (x - V.col(1)).dot(gradients[0]);
};
auto lambda2 = [&V, &gradients](const Eigen::Vector2d& x) -> double {
    return (x - V.col(2)).dot(gradients[1]);
};
auto lambda3 = [&V, &gradients](const Eigen::Vector2d& x) -> double {
    return (x - V.col(0)).dot(gradients[2]);
};

// Store lambda functions in vector with explicit type
std::vector<std::function<double(const Eigen::Vector2d&)>> lambdas = {lambda1, lambda2, lambda3};

Eigen::Matrix3d hard_coded_mat = Eigen::Matrix3d::Zero();

// Use Kahan summation for better numerical precision
for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = i; j < 3; ++j) {
        double sum = 0.0;
        double c = 0.0;  // compensation term
        for(unsigned k = 0; k < quadrature_points.rows(); ++k) {
            Eigen::Vector2d eval_point(
                quadrature_points.row(k).dot(V.row(0)),
                quadrature_points.row(k).dot(V.row(1))
            );
            
            double term = lambdas[i](eval_point) * lambdas[j](eval_point) * 
                         conductivity * quadrature_weights[k] * area;
            
            // Kahan summation algorithm
            double y = term - c;
            double t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        
        hard_coded_mat(i, j) = sum;
        if (i != j) {
            hard_coded_mat(j, i) = sum;
        }
    }
}
std::cout << "hard coded mat : " << std::endl << hard_coded_mat << std::endl ; 


return 0; 
}