#include "DiscreteElasticRods.h"
#include <iostream>

DiscreteElasticRods::DiscreteElasticRods() {}

void DiscreteElasticRods::initSimulation(int nv_, Eigen::VectorXd x_, Eigen::VectorXd theta_, std::vector<bool> is_fixed_, SimParameters params_) {
    nv = nv_;
    x_ref = x_;
    x = x_;
    is_fixed = is_fixed_;
    v.resize(nv*3);
    v.setZero();
    theta = theta_;
//    theta.resize(nv-1);
//    theta.setZero();
    e.resize((nv-1)*3);
    e.setZero();
    length_rest.resize(nv-1);
    length_rest.setZero();
    length.resize(nv-1);
    length.setZero();
    updateEdge();
    updateLength();

    for (size_t i = 0; i < nv-1; i++) {
        length_rest(i) = length(i);
    }
    Eigen::MatrixX3d d3 = unitTangents(x);
    // set up reference frame
    theta.resize(nv-1);
    theta.setZero();
    d1_ref.resize(nv-1,3);
    d2_ref.resize(nv-1,3);

    for (int i = 0; i < nv-1; i++) {
        Eigen::Vector3d d3_i = d3.row(i).transpose();
        d1_ref.row(i) = Eigen::RowVector3d(-d3_i(1), d3_i(0),0);
        Eigen::Vector3d d1_i = d1_ref.row(i).transpose();
        d2_ref.row(i) = (d1_i.cross(d3_i)).transpose();
    }
    d1.resize(nv-1,3);
    d1 = d1_ref;
    d2.resize(nv-1,3);
    d2 = d2_ref;

    kb.resize(nv-2,3);
    updateCurvatureBinormal(d3);
    getMaterialCurvature(kappa_ref);
    twist_rest = getTwist(d2_ref, d3);
    getVoronoiLength(l_ref);

    params = params_;

    vis_gradient.resize(nv);
    vis_stretching_force.resize(nv);
    vis_bending_force.resize(nv);
    vis_twisting_force.resize(nv);
}

void DiscreteElasticRods::simulateOneStep() {
    Eigen::MatrixX3d prev_d3 = unitTangents(x);
    updateCenterlinePosition();
    Eigen::MatrixX3d d3 = unitTangents(x);
    updateMaterialFrame(prev_d3, d3);

    updateEdge();
    updateLength();
    updateCurvatureBinormal(d3);
    Eigen::VectorXd twist = getTwist(d2_ref, d3);

    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> hessian;
    if (verbose) std::cout << "- Computing Energy" << std::endl;
    computeGradientAndHessian(gradient, hessian, d3, twist);
    updateCenterlineVelocity(gradient);
    updateFrameTheta(gradient);
    buildVisualization(gradient);
}

void DiscreteElasticRods::updateCenterlinePosition(void) {
    const double h = params.time_step;
    for (int i = 0; i < nv; i++) {
        if (!is_fixed[i]) x.segment<3>(3*i) += h * v.segment<3>(3*i);
    }
}

void DiscreteElasticRods::updateMaterialFrame(Eigen::MatrixX3d prev_d3, Eigen::MatrixX3d d3) {
    for (int i = 0; i < nv-1; i++) {
        double frame_theta = theta(i);
        d1_ref.row(i) = parallelTransport(d1_ref.row(i), prev_d3.row(i).transpose(), d3.row(i).transpose());
        d2_ref.row(i) = parallelTransport(d2_ref.row(i), prev_d3.row(i).transpose(), d3.row(i).transpose());
        d1.row(i) = std::cos(frame_theta)*d1_ref.row(i) + std::sin(frame_theta)*d2_ref.row(i);
        d2.row(i) = -std::sin(frame_theta)*d1_ref.row(i) + std::cos(frame_theta)*d2_ref.row(i);
    }
}

std::tuple<Eigen::MatrixXd, Eigen::SparseMatrix<double> > DiscreteElasticRods::createZeroGradientAndHessian() {
    int ndof = 3*nv + nv-1;
    Eigen::VectorXd gradient(ndof);
    gradient.setZero();
    Eigen::SparseMatrix<double> hessian(ndof, ndof);
    return std::make_tuple(gradient, hessian);
}

void DiscreteElasticRods::computeGradientAndHessian(Eigen::VectorXd& gradient,
                                                    Eigen::SparseMatrix<double>& hessian,
                                                    Eigen::MatrixX3d& d3,
                                                    Eigen::VectorXd& twist) {
    Eigen::VectorXd stretching_force, bending_force, twisting_force;
    stretching_force.resize(3*nv);
    stretching_force.setZero();
    bending_force.resize(3*nv);
    bending_force.setZero();
    twisting_force.resize(3*nv);
    twisting_force.setZero();

    std::tie(gradient, hessian) = createZeroGradientAndHessian();
    std::vector<Eigen::Triplet<double>> hessian_triplets;

    if (params.stretching_energy_enabled) {
        if (verbose) std::cout << "    - Computing Stretching Energy" << std::endl;
        applyStretchingForce(gradient, hessian_triplets, d3, stretching_force);
    }

    if (params.bending_energy_enabled) {
        if (verbose) std::cout << "    - Computing Bending Energy" << std::endl;
        applyBendingForce(gradient, hessian_triplets, d3, bending_force);
    }

    if (params.twisting_energy_enabled) {
        if (verbose) std::cout << "    - Computing Twisting Energy" << std::endl;
        applyTwistingForce(gradient, hessian_triplets, twist, twisting_force);
    }

    if (params.gravity_enabled) {
        if (verbose) std::cout << "    - Adding Gravity" << std::endl;
        applyGravity(gradient);
    }

    hessian.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());

    buildForceVisualization(stretching_force, bending_force, twisting_force);
}

void DiscreteElasticRods::updateCenterlineVelocity(Eigen::VectorXd& gradient) {
    const double h = params.time_step;
    for (int i = 0; i < nv; i++) {
        // TODO: add mass
        double m = 1;
        v.segment<3>(3*i) -= h * gradient.segment<3>(3 * i) / m;
    }

}

void DiscreteElasticRods::updateFrameTheta(Eigen::VectorXd& gradient) {
    const double h = params.time_step;
    for (int i = 0; i < nv-1; i++) {
        if (!is_fixed[i] || !is_fixed[i+1]) theta(i) -= h * gradient(3*nv+i);
    }
}

void DiscreteElasticRods::updateEdge() {
    for (int i = 0; i < nv-1; i++) {
        e.segment<3>(3*i) = x.segment<3>(3*(i+1))-x.segment<3>(3*i);
    }
}

void DiscreteElasticRods::updateLength() {
    for (int i = 0; i < nv-1; i++) {
        length(i) = e.segment<3>(3*i).norm();
    }
}

void DiscreteElasticRods::updateCurvatureBinormal(Eigen::MatrixX3d d3) {
    for (int i = 0; i < nv-2; i++) {
        Eigen::Vector3d d3_i = d3.row(i).transpose();
        Eigen::Vector3d d3_ip1 = d3.row(i+1).transpose();
        Eigen::Vector3d kb_i = 2*(d3_i).cross(d3_ip1)/(1+d3_i.dot(d3_ip1));
        kb.row(i) = kb_i.transpose();
    }
}

Eigen::MatrixX3d DiscreteElasticRods::unitTangents(const Eigen::VectorXd x_) {
    Eigen::MatrixX3d unit_tangents(nv-1, 3);
    for (int i = 0; i < nv-1; i++)
        unit_tangents.row(i) = ((x_.segment<3>(3*(i+1)) - x_.segment<3>(3*i)).normalized()).transpose();
    return unit_tangents;
}

Eigen::Vector3d DiscreteElasticRods::parallelTransport(Eigen::Vector3d v, Eigen::Vector3d r1, Eigen::Vector3d r2) {
    Eigen::Vector3d k = r1.cross(r2).normalized();
    double theta = std::atan2((r1.cross(r2)).dot(k), r1.dot(r2));
    return v*std::cos(theta) + (k.cross(v))*std::sin(theta) + k*(k.dot(v))*(1-std::cos(theta));
}

double DiscreteElasticRods::applyStretchingForce(Eigen::VectorXd& gradient,
                                                 std::vector<Eigen::Triplet<double>>& hessian,
                                                 Eigen::MatrixX3d& d3,
                                                 Eigen::VectorXd& stretching_force) {
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);
    // ======== compute stretching energy E_s ========
    const double r = params.segment_radius;
    const double k = params.stretching_modulus;
    double E_s = 0.0;
    for (int i = 0; i <= nv-2; i++) {
        const double a_i = r;
        const double b_i = r;
        // TODO: fix A_i = pi * a_j * b_j
        const double A_i = M_PI * a_i * b_i;
        E_s += k*A_i*(std::pow(length(i)-length_rest(i),2)-1)*length(i);
    }
    E_s /= 2;

    // ======== compute stretching force derivative ========
    for (int i = 0; i < nv-1; i++) {
        double e_j_rest = length_rest(i);
        double e_j = length(i);
        Eigen::Vector3d t = d3.row(i).transpose();
        gradient.segment<3>(3 * i) += k * (e_j/e_j_rest-1) * -t * e_j_rest;
        gradient.segment<3>(3 * (i + 1)) += k * (e_j/e_j_rest-1) * t * e_j_rest;
        stretching_force.segment<3>(3 * i) -= k * (e_j/e_j_rest-1) * -t * e_j_rest;
        stretching_force.segment<3>(3 * (i + 1)) -= k * (e_j/e_j_rest-1) * t * e_j_rest;

        // TODO: compute stretching hessian
        // ======== compute bending force hessian ========
    }


    std::vector<Eigen::Triplet<double>> dF_stretching_triplets;

    hessian.insert(hessian.end(), dF_stretching_triplets.begin(), dF_stretching_triplets.end());

    return E_s;
}

void DiscreteElasticRods::getMaterialCurvature(Eigen::MatrixX2d& kappa) {
    kappa.resize(nv-2,2);
    for (int i = 0; i < nv-2; i++) {
        Eigen::VectorXd kb_i = kb.row(i).transpose();
        kappa(i,0) = kb_i.dot((d2.row(i)+d2.row(i+1)).transpose())/2;
        kappa(i,1) = -kb_i.dot((d1.row(i)+d1.row(i+1)).transpose())/2;
    }
}

Eigen::VectorXd DiscreteElasticRods::getTwist(Eigen::MatrixX3d& d2, Eigen::MatrixX3d& d3) {
    Eigen::VectorXd twist(nv-2);
    for (int i = 1; i < nv-1; i++) {
        Eigen::Vector3d vec1 = parallelTransport(d2.row(i-1).transpose(), d3.row(i-1).transpose(), d3.row(i).transpose());
        Eigen::Vector3d vec2 = d2.row(i).transpose();
        Eigen::Vector3d n = (vec1.cross(vec2)).normalized();
        double reference_twist = std::atan2((vec1.cross(vec2)).dot(d3.row(i).transpose()), vec1.dot(vec2));
        std::cout << std::endl;
        twist(i-1) = theta(i) - theta(i-1) + reference_twist;
    }
    return twist;
}

void DiscreteElasticRods::getVoronoiLength(Eigen::VectorXd& l) {
    l.resize(nv-2);
    for (int i = 0; i < nv-2; i++) {
        l(i) = (length_rest(i)+length_rest(i+1))/2;
    }
}

Eigen::Matrix3d DiscreteElasticRods::getCrossMatrix(Eigen::Vector3d v) {
    Eigen::Matrix3d cross_matrix;
    cross_matrix << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;
    return cross_matrix;
}

double DiscreteElasticRods::applyBendingForce(Eigen::VectorXd& gradient,
                                              std::vector<Eigen::Triplet<double>>& hessian,
                                              Eigen::MatrixX3d& d3,
                                              Eigen::VectorXd& bending_force) {
    std::vector<Eigen::Triplet<double>> hess_bending_triplets;
    const double r = params.segment_radius;
    const double E = params.bending_modulus;
    double E_b = 0.0;
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);

    for (int i = 1; i < nv-1; i++) {
        const double a_i = r;
        const double b_i = r;
        // TODO: fix A_i = pi * a_j * b_j
        const double A_i = M_PI * a_i * b_i;
        const double B11 = E * A_i * pow(a_i, 2) / 4;
        const double B22 = E * A_i * pow(b_i, 2) / 4;
        const double l_i = l_ref(i-1);
        const double kappa1_ref_i = kappa_ref(i-1, 0);
        const double kappa2_ref_i = kappa_ref(i-1, 1);
        const double kappa1_i = kappa(i-1, 0);
        const double kappa2_i = kappa(i-1, 1);
        // ======== compute bending energy E_b =========
        E_b += (B11 * std::pow((kappa1_i - kappa1_ref_i), 2) +
                B22 * std::pow((kappa2_i - kappa2_ref_i), 2)) / l_i;

        // ======== compute bending force gradient ========
        Eigen::VectorXd gradient_kappa1_i(11); // 3*3 (x DOF) + 2 (theta DOF)
        Eigen::VectorXd gradient_kappa2_i(11); // 3*3 (x DOF) + 2 (theta DOF)
        const double chi = 1 + d3.row(i-1).dot(d3.row(i));
        Eigen::Vector3d kb_i = kb.row(i-1).transpose();
        Eigen::Vector3d d1_im1 = d1.row(i-1).transpose();
        Eigen::Vector3d d1_i = d1.row(i).transpose();
        Eigen::Vector3d d1_tilde = (d1_i + d1_im1) / chi;
        Eigen::Vector3d d2_im1 = d2.row(i-1).transpose();
        Eigen::Vector3d d2_i = d2.row(i).transpose();
        Eigen::Vector3d d2_tilde = (d2_i + d2_im1) / chi;
        Eigen::Vector3d t_im1 = d3.row(i-1).transpose();
        Eigen::Vector3d t_i = d3.row(i).transpose();
        Eigen::Vector3d t_tilde = (t_i + t_im1) / chi;

        Eigen::Matrix3d de_dx_im1 = -Eigen::Matrix3d::Identity();
        Eigen::Matrix3d de_dx_i = Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dkb_de_im1 = (2 * getCrossMatrix(d3.row(i).transpose()) / chi -
                                      (kb_i / chi) * (d3.row(i-1) + d3.row(i))) / length_rest(i-1);
        Eigen::Matrix3d dkb_de_i = (-2 * getCrossMatrix(d3.row(i-1).transpose()) / chi -
                                    (kb_i / chi) * (d3.row(i-1) + d3.row(i))) / length_rest(i);
        Eigen::Matrix3d test_de_i = (kb_i / chi) * (d3.row(i-1) + d3.row(i)) / length_rest(i);

        Eigen::Matrix3d dkb_dx_im1 = dkb_de_im1 * de_dx_im1;
        Eigen::Matrix3d dkb_dx_i = dkb_de_im1 * de_dx_i + dkb_de_i * de_dx_im1;
        Eigen::Matrix3d dkb_dx_ip1 = dkb_de_i * de_dx_i;
//        gradient_kappa1_i.segment<3>(0) = dkb_dx_im1.transpose() * (d2_im1 + d2_i) / 2;
//        gradient_kappa1_i.segment<3>(3) = dkb_dx_i.transpose() * (d2_im1 + d2_i) / 2;
//        gradient_kappa1_i.segment<3>(6) = dkb_dx_ip1.transpose() * (d2_im1 + d2_i) / 2;
//        gradient_kappa2_i.segment<3>(0) = -dkb_dx_im1.transpose() * (d1_im1 + d1_i) / 2;
//        gradient_kappa2_i.segment<3>(3) = -dkb_dx_i.transpose() * (d1_im1 + d1_i) / 2;
//        gradient_kappa2_i.segment<3>(6) = -dkb_dx_ip1.transpose() * (d1_im1 + d1_i) / 2;
        Eigen::Vector3d dkappa1_i_de_im1 = d2_tilde.cross(-t_i/length_rest(i-1)) - kappa1_i*t_tilde/length_rest(i-1);
        Eigen::Vector3d dkappa1_i_de_i = d2_tilde.cross(t_im1/length_rest(i)) - kappa1_i*t_tilde/length_rest(i);
        Eigen::Vector3d dkappa2_i_de_im1 = d1_tilde.cross(t_i/length_rest(i-1)) - kappa2_i*t_tilde/length_rest(i-1);
        Eigen::Vector3d dkappa2_i_de_i = d1_tilde.cross(-t_im1/length_rest(i)) - kappa2_i*t_tilde/length_rest(i);
        gradient_kappa1_i.segment<3>(0) = dkappa1_i_de_im1.transpose()*de_dx_im1;
        gradient_kappa1_i.segment<3>(3) = dkappa1_i_de_im1.transpose()*de_dx_i + dkappa1_i_de_i.transpose()*de_dx_im1;
        gradient_kappa1_i.segment<3>(6) = dkappa1_i_de_i.transpose()*de_dx_i;
        gradient_kappa2_i.segment<3>(0) = dkappa2_i_de_im1.transpose()*de_dx_im1;
        gradient_kappa2_i.segment<3>(3) = dkappa2_i_de_im1.transpose()*de_dx_i + dkappa2_i_de_i.transpose()*de_dx_im1;
        gradient_kappa2_i.segment<3>(6) = dkappa2_i_de_i.transpose()*de_dx_i;

        gradient_kappa1_i(9) = kb_i.dot(d1_im1) / 2;
        gradient_kappa1_i(10) = kb_i.dot(d1_i) / 2;
        gradient_kappa2_i(9) = -kb_i.dot(d2_im1) / 2;
        gradient_kappa2_i(10) = -kb_i.dot(d2_i) / 2;

        // update dx
        Eigen::Matrix<double, 9, 1> gradient_dx = (B11 * (kappa1_i - kappa1_ref_i) * gradient_kappa1_i.segment<9>(0) + B22 * (kappa2_i - kappa2_ref_i) * gradient_kappa2_i.segment<9>(0)) / l_i;
        gradient.segment<9>(3*(i-1)) += gradient_dx;
        std::cout << "sum: " << gradient_dx.segment<3>(0) + gradient_dx.segment<3>(3) + gradient_dx.segment<3>(6) << std::endl;
        bending_force.segment<9>(3*(i-1)) -= gradient_dx.segment<9>(0);

        // update dtheta
        gradient.segment<2>(3*nv+i-1) += (B11 * (kappa1_i - kappa1_ref_i) * gradient_kappa1_i.segment<2>(9) + B22 * (kappa2_i - kappa2_ref_i) * gradient_kappa2_i.segment<2>(9)) / l_i;

        // TODO: compute bending hessian
        // ======== compute bending force hessian ========
//        Eigen::Matrix<double, 11, 11> hessian_E_b_i;
//        hessian_E_b_i.setZero();
//        Eigen::Matrix<double, 11, 11> hessian_kappa1_i;
//        hessian_kappa1_i.setZero();
//        Eigen::Matrix<double, 11, 11> hessian_kappa2_i;
//        hessian_kappa2_i.setZero();
//
//        for (int j = 3; j < 3; j++) {
//            for (int k = 0; k < j; k++) {
//
//            }
//        }
//
//        hessian_kappa1_i(9, 9) = -kb_i.dot(d2_im1)/2;
//        hessian_kappa1_i(10, 10) = -kb_i.dot(d2_i)/2;
//        hessian_kappa2_i(9, 9) = kb_i.dot(d1_im1)/2;
//        hessian_kappa2_i(10, 10) = kb_i.dot(-d1_i)/2;
//
//        hessian_E_b_i = (B11 * (gradient_kappa1_i * gradient_kappa1_i.transpose() + (kappa1_i - kappa1_ref_i) * hessian_kappa1_i)
//                        +B22 * (gradient_kappa2_i * gradient_kappa2_i.transpose() + (kappa2_i - kappa2_ref_i) * hessian_kappa2_i))/l_i;
//
//        // update dx_a dx_b hessian
//        for (int j = 0; j < 3; j++) {
//            hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*(i-1)+j, 3*(i-1), hessian_E_b_i(j, 0)));
//            hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*(i-1)+j, 3*(i-1)+1, hessian_E_b_i(j, 1)));
//            hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*(i-1)+j, 3*(i-1)+2, hessian_E_b_i(j, 2)));
//        }
//
//        // update dx_a dtheta_a & dtheta_a dx_a hessian
//        for (int j = 0; j < 2; j++) {
//            for (int k = 0; k < 3; k++) {
//                double hessian_entry = hessian_E_b_i(3+j, k);
//                hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*nv+i-1+j, 3*(i-1)+k, hessian_entry));
//                hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*(i-1)+k, 3*nv+i-1+j, hessian_entry));
//            }
//        }
//
//        // update dtheta_a dtheta_b hessian
//        for (int j = 0; j < 2; j++) {
//            hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*nv+i-1+j, 3*nv+i-1, hessian_E_b_i(3+j, 3+j)));
//            hess_bending_triplets.emplace_back(Eigen::Triplet<double>(3*nv+i-1+j, 3*nv+i, hessian_E_b_i(3+j, 3+j+1)));
//        }
    }
    E_b /= 2;
    hessian.insert(hessian.end(), hess_bending_triplets.begin(), hess_bending_triplets.end());
    return E_b;
}

double DiscreteElasticRods::applyTwistingForce(Eigen::VectorXd& gradient,
                                               std::vector<Eigen::Triplet<double>>& hessian,
                                               Eigen::VectorXd& twist,
                                               Eigen::VectorXd& twisting_force) {
    std::vector<Eigen::Triplet<double>> hess_twisting_triplets;
    const double r = params.segment_radius;
    const double G = params.twisting_modulus;
    double E_t = 0.0;
    Eigen::MatrixX2d kappa;
    getMaterialCurvature(kappa);
    for (int i = 1; i < nv-1; i++) {
        const double a_i = r;
        const double b_i = r;
        //TODO: fix A_i = pi * a_j * b_j
        const double A_i = M_PI * a_i * b_i;
        double beta_i = G * A_i * (std::pow(a_i, 2) + std::pow(b_i, 2)) / 4;
        const double l_i = l_ref(i - 1);
        // ======== compute twisting energy E_t ========
        E_t += beta_i / l_i * std::pow((twist(i-1) - twist_rest(i-1)), 2);

        // ======== compute twisting force gradient ========
        Eigen::VectorXd gradient_m_i(11); // 3*3 (x DOF) + 2 (theta DOF)

        Eigen::Vector3d kb_i = kb.row(i-1).transpose();
        double e_im1 = length(i-1);
        double e_i = length(i);
        // update dx
        gradient_m_i.segment<3>(0) = -kb_i / (2*e_im1);
        gradient_m_i.segment<3>(3) = -kb_i / (2*e_i) + kb_i / (2*e_im1);
        gradient_m_i.segment<3>(6) = kb_i / (2*e_i);
        // update dtheta
        gradient_m_i(9) = -1.;
        gradient_m_i(10) = 1.;

        Eigen::Matrix<double, 9, 1> gradient_dx = beta_i / l_i * (twist(i-1) - twist_rest(i - 1)) * gradient_m_i.segment<9>(0);
        gradient.segment<9>(3*(i-1)) += gradient_dx;
        twisting_force.segment<9>(3*(i-1)) -= gradient_dx;
        gradient.segment<2>(3*nv+i-1) += beta_i / l_i * (twist(i-1) - twist_rest(i - 1)) * gradient_m_i.segment<2>(9);

        // TODO: compute twisting Hessian dF
        // ======== compute bending force hessian ========
    }
    E_t /= 2;
    hessian.insert(hessian.end(), hess_twisting_triplets.begin(), hess_twisting_triplets.end());
    return E_t;
}

double DiscreteElasticRods::applyGravity(Eigen::VectorXd& gradient) {
    const double g = params.gravity_G;
    for (int i = 0; i < nv; i++) {
        // TODO: add mass
        double m = 1;
        gradient(3 * i + 1) -= m * g;
    }
    return 0;
}

void DiscreteElasticRods::buildVisualization(Eigen::VectorXd& gradient) {
    for (int i = 0; i < nv; i++) {
        vis_gradient[i] = gradient.segment<3>(3*i);
    }
}

void DiscreteElasticRods::buildForceVisualization(Eigen::VectorXd& stretching_force,
                                                  Eigen::VectorXd& bending_force,
                                                  Eigen::VectorXd& twisting_force) {
    for (int i = 0; i < nv; i++) {
        vis_stretching_force[i] = stretching_force.segment<3>(3*i);
        vis_bending_force[i] = bending_force.segment<3>(3*i);
        vis_twisting_force[i] = twisting_force.segment<3>(3*i);
    }
}
