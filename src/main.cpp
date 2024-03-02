#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "eigenmvn.h"
#include "matern.h"
#include "calc_posterior/posterior.h"
#include "chrono"
#include "coordinates/coordinates.h"
#include "ar_model/ar_class.h"
#include<random>
#include"debug_functions/debug.h"
// #include"build/proto/ydata.pb.h"
#include"build/proto/parsedata.pb.h"
#include"protocpp/serialize.h"



// Source files

void print_vec(std::vector<double>& vector);

int main(int argc,char* argv[]) {

    // argv[1] is the path to parsed data file

    //argv[2] is nu
    double nu = std::stod(argv[2]);
    //argv[3] & argv[4] are prior a_phi and b_phi from user input
    double first = std::stod(argv[3]);
    double second = std::stod(argv[4]);
    std::pair<double, double> ab_phi_prior = {first, second};
    //argv[5] is n_iter
    unsigned int n_iter = std::stoi(argv[5]);
    //argv[6] is burn_in
    unsigned int burn_in = std::stoi(argv[6]);

    // algo options
    u_int64_t seed = 112083918;
    bool b = true;

    ///////// DATA PARSING ///////////
    std::cout << "Parsing data..." << std::endl;
    // // // Deserialize the data
    parsedata::input_data data;
    
    std::fstream input(argv[1], std::ios::in | std::ios::binary);
    if (!data.ParseFromIstream(&input)) {
        data.PrintDebugString();
        std::cerr << "Failed to parse data." << std::endl;
        return 1;
        }
    else{
        std::cout << "Data was parsed successfully. Continuing ..." << std::endl;
    }
    // Accessing the deserialized data
    parsedata::vector y_parse = data.y();
    parsedata::matrix x_parse = data.x();
    google::protobuf::RepeatedPtrField<parsedata::location> loc = data.loc();


    // Convert protobuf data structures to C++ data structures
    std::vector<std::vector<double>> xdata;
    for (const parsedata::vector& col : x_parse.m_vec()) {
        std::vector<double> colValues(col.vec_value().begin(), col.vec_value().end());
        xdata.push_back(colValues);
    }

    std::vector<double> ydata;
    for (const double& vec_value : y_parse.vec_value()) {
            ydata.push_back(vec_value);

    }

    const unsigned int p = xdata.size();
    std::vector<coord> coord_store_vec;

    for (const parsedata::location& locEntry : loc) {
        coord temp(locEntry.lat(), locEntry.long_());
        coord_store_vec.push_back(temp);
    }

    const unsigned int N = coord_store_vec.size();
    const unsigned int T = ydata.size()/N;
    // divide y
    std::vector<Eigen::VectorXd> Y;
    for (auto it = ydata.begin(); it != ydata.end(); it += N) {
        std::vector<double> y_sliced(it, it + N);
        Eigen::Map<Eigen::VectorXd> temp(y_sliced.data(), N);
        Y.push_back(temp);
    }

    // divide X
    std::vector<Eigen::MatrixXd> X;

        for (int i = 0; i < N*T; i += N) {
            std::vector<double> temp;
            temp.reserve(N*T);
            for (int j = 0; j < p; ++j) {
                temp.insert(temp.end(), xdata[j].begin()+ i, xdata[j].begin()+N + i);
            }
            Eigen::Map<Eigen::MatrixXd> x_t(temp.data(), N, p);
            X.push_back(x_t);
        }

    ///////// DATA PARSING END ////////

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ar_model a(Y, X, coord_store_vec, ab_phi_prior, nu);
    a.init();
    for(int i = 0; i <burn_in; ++i) {
        a.sample();
        if(i % 100 == 0) {
            std::cout<< "BURN-IN: Iteration " << i << " finished" << std::endl;
       }

    }
    std::cout << "Burn in finished. Continuing..." << std::endl;
    for (int i = 0; i < n_iter; ++i) {
        a.sample();
        a.write_curr_state();
        a.track_pmcc();
       if(i % 100 == 0) {
            std::cout << "Iteration " << i << " finished" << std::endl;
       }

    }
    std::cout << "Sampling finished. Importing data from C++..." << std::endl;
    a.serialize();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    std::cout << "phi acceptance rate: " << a.get_acceptance_rate() << std::endl;
    std::cout << "PMCC: " <<a.calc_pmcc() << std::endl;
    return 0;
}
