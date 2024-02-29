//
// Created by Daniel Adamovic on 21/08/23.
//

#ifndef AR_GIBBS_COORDINATES_CPP
#define AR_GIBBS_COORDINATES_CPP

#include<iostream>
#include <cmath>
#include<Eigen/Dense>
#include<../matern.h>

class coord {
private:
    double x;
    double y;
public:
    coord() = default;
    coord(double x, double y): x(x), y(y){};
    double get_x() const {
        return x;
    };
    double get_y() const{
        return y;
    };
    void set_x(double& x){
        this->x = x;
    };
    void set_y(double& y) {
        this->y = y;
    }
};

double eucl_dist(coord& c1, coord& c2);

Eigen::MatrixXd eucl_dist_matrix(std::vector<coord>& coordinates);

#endif //AR_GIBBS_COORDINATES_CPP
