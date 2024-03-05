#ifndef MATERN_H
#define MATERN_H

#include <iostream>
#include<cmath>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include "Eigen/Dense"
#include <numbers>


// distance function and Gamma function definition
// bessel function from boost

static double matern(double& x, double& phi, double& nu){
    double ans;
    if(x == 0)
       ans = 1;
    else
       ans = (1.0 / ( exp(std::lgamma(nu)) * pow(2, nu - 1.0) )) * pow(sqrt(2*nu)*x*phi, nu) * boost::math::cyl_bessel_k(nu, sqrt(2*nu)*x*phi);
    return ans;
 }

static Eigen::MatrixXd calc_matern_mat(Eigen::MatrixXd& coord_mat, double& phi, double& nu){
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(coord_mat.cols(), coord_mat.cols());
    for (int i = 0; i < coord_mat.cols(); ++i) {
        for (int j = 0; j < coord_mat.cols(); ++j) {
            temp(i,j) = matern(coord_mat(i,j), phi, nu);
        }
    }
    return temp;
}
#endif

