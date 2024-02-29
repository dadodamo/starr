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
    std::cout << argv[1] << std::endl;

    //argv[2] is nu
    // std::cout << argv[2] << std::endl;

    //argv[3] & argv[4] are prior a_phi and b_phi from user input
    // std::cout <<"a_phi is: " <<argv[3] << " and b_phi: " << argv[4]<< std::endl;

    // algo options
    u_int64_t seed = 112083918;
    bool b = true;
    double nu = 0.5;

    ///////// DATA PARSING ///////////

    // std::ifstream input(argv[1], std::ios::binary);
    // if (!input) {
    //     std::cerr << "Failed to open dataparsed.bin" << std::endl;
    //     return 1;
    // }
    // std::stringstream buffer;
    // buffer << input.rdbuf();
    // std::string serialized_data = buffer.str();

    // // Deserialize the data
    parsedata::input_data data;
    // if (!data.ParseFromString(serialized_data)) {
    //     std::cerr << "Failed to parse serialized data" << std::endl;
    //     return 1;
    // }
    std::fstream input(argv[1], std::ios::in | std::ios::binary);
    if (!data.ParseFromIstream(&input)) {
        data.PrintDebugString();
        std::cerr << "Failed to parse address book." << std::endl;
        return 1;
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
    std::cout << N << std::endl;
    std::cout << T << std::endl;
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
//    Eigen::Map<Eigen::VectorXd> eigenVector(stdVector.data(), stdVector.size());

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ar_model a(Y, X, coord_store_vec, nu);
    a.init();
    a.standardize();
    unsigned int n_iter = 5000;
    unsigned int burn_in = 1000;
    for(int i = 0; i <burn_in; ++i) {
        a.sample();
        if(i % 100 == 0) {
            std::cout<< "BURN-IN: Iteration " << i << " finished" << std::endl;
       }

    }
    for (int i = 0; i < n_iter; ++i) {
        a.sample();
        a.write_curr_state();
        a.track_pmcc();
       if(i % 100 == 0) {
            std::cout << "Iteration " << i << " finished" << std::endl;
       }

    }
    std::cout << "acceptance rate: " << a.get_acceptance_rate() << std::endl;
    std::cout << a.calc_pmcc() << std::endl;
    a.serialize();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    return 0;
}
