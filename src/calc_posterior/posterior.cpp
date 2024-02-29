#include "posterior.h"
Eigen::VectorXd post::calc_mean_beta( const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec, double& rho){
        Eigen::VectorXd mean(x_store_vec[0].cols());
        for (int t = 1; t < x_store_vec.size(); t++) {
            mean += x_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - rho*ot_store_vec[t-1]);
        }
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_beta(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, const double& sigma_beta_prior){
        auto id_mat = Eigen::MatrixXd::Identity(x_store_vec[0].cols(),x_store_vec[0].cols() );
        Eigen::MatrixXd temp(x_store_vec[0].cols(), x_store_vec[0].cols());
        //iteration starts at first element of x_store_vec which is t = 1
        for (int t = 0; t < x_store_vec.size()-1; ++t) {
            temp += x_store_vec[t].transpose() * covar_w_inv * x_store_vec[t];
        }
        temp +=  id_mat * 1/sigma_beta_prior;
        return temp.inverse();
    };

// rho
    double post::calc_mean_rho(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::VectorXd& beta){
        double mean = 0;
        for (int t = 1; t < ot_store_vec.size(); ++t) {
            mean += ot_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - x_store_vec[t-1]*beta);
        }
        return mean;
    };
    double post::calc_var_rho(std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::MatrixXd& covar_w_inv,const double& sigma_rho_prior){
        double var = 0.;
        for (int t = 1; t < ot_store_vec.size(); ++t) {
            var += ot_store_vec[t-1].transpose() * covar_w_inv * ot_store_vec[t-1];
        }
        return 1./(var + 1./sigma_rho_prior);
    };

// O_T
    Eigen::VectorXd post::calc_mean_eff_T( const Eigen::VectorXd& y_T, const Eigen::MatrixXd& x_T, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev,
                                           Eigen::VectorXd& beta, double& rho, double& sigma_eps){
        Eigen::VectorXd mean(y_T.rows());
        mean = (y_T * 1/sigma_eps)  + covar_w_inv*(rho*o_prev + x_T*beta);
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_eff_T(double& sigma_eps, Eigen::MatrixXd& covar_w_inv){
        auto id_mat = Eigen::MatrixXd::Identity(covar_w_inv.cols(),covar_w_inv.cols() );
        Eigen::MatrixXd cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (id_mat * 1/sigma_eps + covar_w_inv);
        return cov.inverse();
    };

// O_0
    Eigen::VectorXd post::calc_mean_eff_0(const Eigen::MatrixXd&  x_1, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_1, Eigen::VectorXd& beta,const Eigen::MatrixXd& S_0_inv,
                                    Eigen::VectorXd& mu_0, double& rho, double& sigma_0){
        Eigen::VectorXd mean =  ( S_0_inv * mu_0 * 1/sigma_0) + covar_w_inv * rho* (o_1 - x_1 * beta);
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_eff_0(Eigen::MatrixXd& covar_w_inv,const Eigen::MatrixXd& S_0_inv, double& rho, double& sigma_0){
        Eigen::MatrixXd cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = ( (pow(rho, 2) *  covar_w_inv) + (1./sigma_0 * S_0_inv));
        return cov.inverse();
    };

// O_t
    Eigen::VectorXd post::calc_mean_eff_t(const Eigen::VectorXd& y_curr,const Eigen::MatrixXd& x_curr,const Eigen::MatrixXd& x_next, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev, Eigen::VectorXd& o_next,
                                           Eigen::VectorXd& beta, double& rho, double& sigma_eps) {
        Eigen::VectorXd mean = y_curr/sigma_eps + covar_w_inv * (x_curr*beta + rho  * (o_prev + o_next - x_next*beta));
        return mean;
    };
// note: It is the exact same as for T, hence I will probably just use one declaration of the function
    Eigen::MatrixXd post::calc_cov_eff_t(double& sigma_eps, Eigen::MatrixXd& covar_w_inv, double& rho){
        auto id_mat = Eigen::MatrixXd::Identity(covar_w_inv.cols(),covar_w_inv.cols());
        Eigen::MatrixXd cov = (id_mat/sigma_eps + (1+pow(rho,2)) * covar_w_inv);
        return cov.inverse();
    };

    std::pair<double, double> post::calc_a_b_sigma_eps(const double& a_prior,const double& b_prior,
                                               const unsigned int& n,const  unsigned int& T ,const std::vector<Eigen::VectorXd>& y_store_vec, std::vector<Eigen::VectorXd>& o_store_vec){
        auto id_mat = Eigen::MatrixXd::Identity(y_store_vec[0].rows(), y_store_vec[0].rows());
        std::pair<double, double> param;
        param.first = a_prior + 0.5*n*T ;
        double temp = 0;
        for (int t = 1; t < o_store_vec.size(); t++) {
            temp += (y_store_vec[t-1] - o_store_vec[t]).transpose() * (y_store_vec[t-1] - o_store_vec[t]);
        }
        param.second = b_prior +  0.5*temp;
        return param;
    };
    std::pair<double, double> post::calc_a_b_sigma_w(const double& a_prior,const double& b_prior,
                                             const unsigned int& n,const unsigned  int& T ,const std::vector<Eigen::MatrixXd>& x_store_vec, const Eigen::MatrixXd& matern_inv ,
                                             std::vector<Eigen::VectorXd>& o_store_vec, Eigen::VectorXd& beta, double& rho){
        std::pair<double, double> param;
        param.first = a_prior + 0.5 * n*T;
        double temp = 0;
        for (int t = 1; t <o_store_vec.size(); t++) {

            temp += (o_store_vec[t] - rho *o_store_vec[t-1]  - x_store_vec[t-1] * beta).transpose() * matern_inv *(o_store_vec[t] - rho *o_store_vec[t-1]  - x_store_vec[t-1] * beta);
        }
        param.second = b_prior + 0.5* temp;
        return param;
    };
    std::pair<double, double> post::calc_a_b_sigma_0(const double& a_prior, const double& b_prior,
                                             const unsigned int& n, const Eigen::MatrixXd& S_0_inv ,Eigen::VectorXd& o_0,Eigen::VectorXd& mu_0){
        std::pair<double, double> param;
        param.first = a_prior + 0.5*n;
        double temp = (o_0 - mu_0).transpose() * S_0_inv * (o_0 - mu_0);
        param.second = b_prior + 0.5*temp;
        return param;
    };
    Eigen::VectorXd post::calc_mean_mu_0(const Eigen::MatrixXd& S_0_inv, Eigen::VectorXd& o_0, double& sigma_0){
        return (S_0_inv * o_0)/sigma_0;
    };
    Eigen::MatrixXd post::calc_cov_mu_0(const Eigen::MatrixXd& S_0_inv, double& sigma_0, const double& sigma_mu_prior){
        auto id_mat = Eigen::MatrixXd::Identity(S_0_inv.cols(),S_0_inv.cols());
        Eigen::MatrixXd cov =  S_0_inv/sigma_0 +  id_mat/sigma_mu_prior;
        return cov.inverse();
    };
    double post::target_ratio_phi(double& curr_phi, double& cand_phi, Eigen::MatrixXd& coord_mat, std::vector<Eigen::VectorXd> o_store,const std::vector<Eigen::MatrixXd> x_store, Eigen::VectorXd& beta,
                                  Eigen::VectorXd& mu_0, double& rho, double& sigma_w, double& sigma_0, std::pair<double, double>& ab_phi_prior, double& nu) {
        Eigen::MatrixXd cand_matern_mat = calc_matern_mat(coord_mat, cand_phi, nu);
        Eigen::MatrixXd cand_matern_inv = cand_matern_mat.inverse();
        Eigen::MatrixXd curr_matern_mat = calc_matern_mat(coord_mat, curr_phi, nu);
        Eigen::MatrixXd curr_matern_inv = curr_matern_mat.inverse();
        double sum_cand = 0;
        double sum_curr = 0;
        for (int t = 1; t < x_store.size(); ++t) {


            sum_cand +=  (o_store[t] - rho*o_store[t-1] - x_store[t-1]*beta).transpose()*cand_matern_inv *(o_store[t] - rho*o_store[t-1] - x_store[t-1]*beta);
            sum_curr += (o_store[t] - rho*o_store[t-1] - x_store[t-1]*beta).transpose()*curr_matern_inv *(o_store[t] - rho*o_store[t-1] - x_store[t-1]*beta);
        }

        double log_cand_lkhd = (ab_phi_prior.first - 1)*log(cand_phi) - ab_phi_prior.second*cand_phi  - (x_store.size()/2.)*log(cand_matern_mat.determinant()) - (1/(2*sigma_w)) * sum_cand -0.5* log(cand_matern_mat.determinant())
                 -1/(2*sigma_0) *(o_store[0] - mu_0).transpose() *cand_matern_inv*(o_store[0] - mu_0);

//        double log_cand_lkhd =  - (x_store.size()/2.)*log(cand_matern_mat.determinant()) - (1/(2*sigma_w)) * sum_cand -0.5* log(cand_matern_mat.determinant())
//                               -1/(2*sigma_0) *(o_store[0] - mu_0).transpose() *cand_matern_inv*(o_store[0] - mu_0);

        double log_curr_lkhd =  (ab_phi_prior.first - 1)*log(curr_phi) - ab_phi_prior.second*curr_phi - (x_store.size()/2.)*log(curr_matern_mat.determinant())  -(1/(2*sigma_w)) * sum_curr -0.5* log(curr_matern_mat.determinant())
                -1/(2*sigma_0) *(o_store[0] - mu_0).transpose() *curr_matern_inv*(o_store[0] - mu_0);
//
//        double log_curr_lkhd =   - (x_store.size()/2.)*log(curr_matern_mat.determinant())  - (1/(2*sigma_w)) * sum_curr -0.5* log(curr_matern_mat.determinant())
//                                -(1/(2*sigma_0)) *(o_store[0] - mu_0).transpose() *curr_matern_inv*(o_store[0] - mu_0);

        return exp(log_cand_lkhd - log_curr_lkhd);
    }

