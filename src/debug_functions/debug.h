//
// Created by Daniel Adamovic on 17/10/23.
//

#include <iostream>
#include<cmath>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include "Eigen/Dense"

#ifndef AR_GIBBS_DEBUG_H
#define AR_GIBBS_DEBUG_H

static void check_eigenvalues(Eigen::MatrixXd mat) {
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(mat);
    std::complex<double> E;
    for(int i = 0; i < mat.rows(); ++i) {
        E = eigensolver.eigenvalues().col(0)[i];
        std::cout << E << std::endl;
    }
}

#endif //AR_GIBBS_DEBUG_H
