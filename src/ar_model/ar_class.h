#ifndef AR_GIBBS_AR_CLASS_H
#define AR_GIBBS_AR_CLASS_H

#include <iostream>
#include "Eigen/Dense"
#include "eigenmvn.h"
#include "../coordinates/coordinates.h"
#include "matern.h"
#include "../build/proto/ydata.pb.h"
#include "../build/proto/paramdata.pb.h"
#include "../protocpp/serialize.h"

#include <fstream>
#include<string>
#include "../calc_posterior/posterior.h"
#include <random>
#include "../debug_functions/debug.h"

// DEBUG MODE


class ar_model {
private:
    //dimensions
    unsigned int T;
    unsigned int N;
    unsigned int p;

    // data passed to class
    std::vector<Eigen::VectorXd> y;
    std::vector<Eigen::MatrixXd> X;
    std::vector<coord> coordinates;
    double nu;

    //sampled parameters
    std::vector<Eigen::VectorXd> o_store;
    Eigen::VectorXd beta;
    double rho;
    Eigen::VectorXd mu_0;
    double phi;
    double sigma_eps;
    double sigma_w;
    double sigma_0;


    // priors chosen
    double beta_sig_prior = 1;
    double rho_sig_prior = 1.;
    double mu0_sig_prior = 1.;

    //inverse gamma group
    std::pair<double, double> ab_eps_prior = {2,4};
    std::pair<double, double> ab_w_prior = {2,4};
    std::pair<double, double> ab_0_prior = {2,4};
    // phi prior and candidate variance
    std::pair<double, double> ab_phi_prior;
    double phi_cand_var = 1;

    // matrices
    Eigen::MatrixXd coord_mat;
    Eigen::MatrixXd matern_cov;
    Eigen::MatrixXd matern_inv;
    Eigen::MatrixXd w_full_cov_inv;

    // algo options
    bool use_cholesky = true;
    u_int64_t seed = 1;
    std::mt19937 generator = std::mt19937(1);

    // phi MH step variables
    double phi_accept_count = 0;
    unsigned int batch_count_50 = 1;
    double phi_s = 0;

    //iter count
    unsigned int iter_count = 0;
    unsigned int inclburn_iter_count = 0;
    unsigned int phi_accept_batch_count = 0;

    // PMCC variables
    std::normal_distribution<double> pmcc_y_sampler = std::normal_distribution<double>(0,1);
    Eigen::VectorXd sampled_y;
    Eigen::VectorXd sampled_y_sum;
    Eigen::VectorXd sampled_y_sum_sq;
    Eigen::VectorXd fitted_y_values;


    //samplers
    Eigen::EigenMultivariateNormal<double> o_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1,1), false, 1);
    Eigen::EigenMultivariateNormal<double> beta_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1),Eigen::MatrixXd::Identity(1,1), false, 1);
    Eigen::EigenMultivariateNormal<double> mu0_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1,1), false, 1);
    std::normal_distribution<double> rho_sampler = std::normal_distribution<double>(0,1);
    std::gamma_distribution<double> sig_eps_sampler = std::gamma_distribution<double>(1, 1);
    std::gamma_distribution<double> sig_w_sampler = std::gamma_distribution<double>(1,1);
    std::gamma_distribution<double> sig_0_sampler = std::gamma_distribution<double>(1,1);
    std::normal_distribution<double> phi_sampler = std::normal_distribution<double>(0,1);
    std::uniform_real_distribution<double> unif = std::uniform_real_distribution<double>(0,1);

    sampler_data::samples sample_stream;
    y_data::full_y y_stream;

    // NA handling
    std::vector<std::pair<int, int>> na_values;

public:
    ar_model(
            std::vector<Eigen::VectorXd>& y_store,
            std::vector<Eigen::MatrixXd>& x_store,
            std::vector<coord>& coord_vec,
            std::pair<double, double>& ab_phi_prior,
            double& nu
    );

    void init();
    void sample();
    void write_curr_state();

    void serialize();
    void track_pmcc();
    double calc_pmcc();

    double get_acceptance_rate(){
        return phi_accept_count/inclburn_iter_count;
    }
};

#endif
