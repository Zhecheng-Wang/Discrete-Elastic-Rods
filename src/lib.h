#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct context {
    std::vector<Eigen::Matrix<double, 2, 2>> B;
    Eigen::VectorXd A;
    double E_b;
    Eigen::VectorXd kappa_1;
    Eigen::VectorXd kappa_2;
    Eigen::VectorXd barkappa_1;
    Eigen::VectorXd barkappa_2;
    double E_t;
    Eigen::VectorXd beta;
    double E_s;

    context(
        const std::vector<double> & a,
        const std::vector<double> & b,
        const std::vector<double> & barl,
        const double & E,
        const std::vector<Eigen::Matrix<double, 3, 1>> & kappab,
        const std::vector<Eigen::Matrix<double, 3, 1>> & barkappab,
        const std::vector<Eigen::Matrix<double, 3, 1>> & d_1,
        const std::vector<Eigen::Matrix<double, 3, 1>> & d_2,
        const std::vector<Eigen::Matrix<double, 3, 1>> & tilded_1,
        const std::vector<Eigen::Matrix<double, 3, 1>> & tilded_2,
        const std::vector<Eigen::Matrix<double, 3, 1>> & bard_1,
        const std::vector<Eigen::Matrix<double, 3, 1>> & bard_2,
        const std::vector<Eigen::Matrix<double, 3, 1>> & bartilded_1,
        const std::vector<Eigen::Matrix<double, 3, 1>> & bartilded_2,
        const double & G,
        const std::vector<double> & m,
        const std::vector<double> & barm,
        const std::vector<Eigen::Matrix<double, 3, 1>> & e,
        const std::vector<Eigen::Matrix<double, 3, 1>> & bare,
        const double & k_s)
    {
        const long dim_0 = a.size();
        const long dim_1 = e.size();
        assert( b.size() == dim_0 );
        assert( barl.size() == dim_0 );
        assert( kappab.size() == dim_0 );
        assert( barkappab.size() == dim_0 );
        assert( d_1.size() == dim_0 );
        assert( d_2.size() == dim_0 );
        assert( tilded_1.size() == dim_0 );
        assert( tilded_2.size() == dim_0 );
        assert( bard_1.size() == dim_0 );
        assert( bard_2.size() == dim_0 );
        assert( bartilded_1.size() == dim_0 );
        assert( bartilded_2.size() == dim_0 );
        assert( m.size() == dim_0 );
        assert( barm.size() == dim_0 );
        assert( bare.size() == dim_1 );
        // A_i = π a_i b_i
        A.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            A[i-1] = M_PI * a.at(i-1) * b.at(i-1);
        }
        // B_i = EA_i/4 [a_i^2 0; 0 b_i^2]
        B.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            Eigen::Matrix<double, 2, 2> B_i_0;
            B_i_0 << pow(a.at(i-1), 2), 0,
            0, pow(b.at(i-1), 2);
            B.at(i-1) = E * A[i-1] / double(4) * B_i_0;
        }
        // `kappa_1`_i = `kappab`_i ⋅ (`\tilde{d_2}`_i + `d_2`_i)/2
        kappa_1.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            kappa_1[i-1] = (kappab.at(i-1)).dot((tilded_2.at(i-1) + d_2.at(i-1))) / double(2);
        }
        // `kappa_2`_i = -`kappab`_i ⋅ (`\tilde{d_1}`_i + `d_1`_i)/2
        kappa_2.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            kappa_2[i-1] = -(kappab.at(i-1)).dot((tilded_1.at(i-1) + d_1.at(i-1))) / double(2);
        }
        // `\bar{kappa}_1`_i = `\bar{kappab}`_i ⋅ (`\bar{\tilde{d_2}}`_i + `\bar{d_2}`_i)/2
        barkappa_1.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            barkappa_1[i-1] = (barkappab.at(i-1)).dot((bartilded_2.at(i-1) + bard_2.at(i-1))) / double(2);
        }
        // `\bar{kappa}_2`_i = -`\bar{kappab}`_i ⋅ (`\bar{\tilde{d_1}}`_i + `\bar{d_1}`_i)/2
        barkappa_2.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            barkappa_2[i-1] = -(barkappab.at(i-1)).dot((bartilded_1.at(i-1) + bard_1.at(i-1))) / double(2);
        }
        double sum_4 = 0;
        for(int i=1; i<=barl.size(); i++){
            sum_4 += 1 / double(barl.at(i-1)) * (B.at(i-1)(1-1, 1-1) * pow((kappa_2[i-1] - barkappa_2[i-1]), 2) + B.at(i-1)(2-1, 2-1) * pow((kappa_1[i-1] - barkappa_1[i-1]), 2));
        }
        // `E_b` = 1/2 ∑_i 1/`\bar{l}`_i(B_i,1,1 (`kappa_2`_i - `\bar{kappa}_2`_i)^2 + B_i,2,2 (`kappa_1`_i - `\bar{kappa}_1`_i)^2)
        E_b = 1 / double(2) * sum_4;
        // beta_i = GA_i(a_i^2+b_i^2)/4
        beta.resize(dim_0);
        for( int i=1; i<=dim_0; i++){
            beta[i-1] = G * A[i-1] * (pow(a.at(i-1), 2) + pow(b.at(i-1), 2)) / double(4);
        }
        double sum_6 = 0;
        for(int i=1; i<=beta.size(); i++){
            sum_6 += beta[i-1] / double(barl.at(i-1)) * pow((m.at(i-1) - barm.at(i-1)), 2);
        }
        // `E_t` = 1/2 ∑_i beta_i/`\bar{l}`_i(m_i - `\bar{m}`_i)^2
        E_t = 1 / double(2) * sum_6;
        double sum_7 = 0;
        for(int j=1; j<=e.size(); j++){
            sum_7 += k_s * pow(((e.at(j-1)).lpNorm<2>() / double((bare.at(j-1)).lpNorm<2>()) - 1), 2) * (bare.at(j-1)).lpNorm<2>();
        }
        // `E_s` = 1/2 ∑_j `k_s` (‖e_j‖ / ‖`\bar{e}`_j‖ - 1)^2 ‖`\bar{e}`_j‖
        E_s = 1 / double(2) * sum_7;
    }
};

