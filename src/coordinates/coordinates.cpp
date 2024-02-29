//
// Created by Daniel Adamovic on 26/09/23.
//
#include "coordinates.h"

double eucl_dist(coord& c1, coord& c2) {
    return sqrt(pow(c1.get_x()-c2.get_x(),2) + pow(c1.get_y() - c2.get_y(), 2));
}

Eigen::MatrixXd eucl_dist_matrix(std::vector<coord>& coordinates){
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(coordinates.size(), coordinates.size());
    for (int i = 0; i < coordinates.size(); ++i) {
        for (int j = 0; j < coordinates.size(); ++j) {
            double dist = eucl_dist(coordinates[i], coordinates[j]) ;
            temp(i, j) = dist;
        }
    }
    return temp;
}