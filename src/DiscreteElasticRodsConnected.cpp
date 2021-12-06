#include "DiscreteElasticRodsConnected.h"
#include <iostream>

void Connector::initConnector(std::vector<std::shared_ptr<DERC> > rods,
        int _rod_index,
        std::pair<int, double>& node_index_segment_ratio_pair,
        int _other_rod_index,
        std::pair<int, double>& other_node_index_segment_ratio_pair,
        double _k)
{
    rod_index = _rod_index;
    std::tie(node_index, segment_ratio) = node_index_segment_ratio_pair;
    other_rod_index = _other_rod_index;
    std::tie(other_node_index, other_segment_ratio) = other_node_index_segment_ratio_pair;
    k = _k;
    Eigen::Vector3d pos = rods[rod_index]->x.segment<3>(3*node_index);
    Eigen::VectorXd& other_rod_nodes = rods[other_rod_index]->x;
    Eigen::Vector3d orig = other_rod_nodes.segment<3>(3*other_rod_index);
    Eigen::Vector3d other_pos = orig + (orig-other_rod_nodes.segment<3>(3*other_rod_index-1)).normalized()*other_segment_ratio;
    length = rest_length = (pos-other_pos).norm();
}

void Connector::computeConnectorEnergy(std::vector<std::shared_ptr<DERC> > rods)
{
    Eigen::Vector3d pos = rods[rod_index]->x.segment<3>(3*node_index);
    Eigen::VectorXd& other_rod_nodes = rods[other_rod_index]->x;
    Eigen::Vector3d orig = other_rod_nodes.segment<3>(3*other_rod_index);
    Eigen::Vector3d other_pos = orig + (orig-other_rod_nodes.segment<3>(3*other_rod_index-1)).normalized()*other_segment_ratio;
    length = (pos-other_pos).norm();
    springEnergy = k/2*std::pow(length-rest_length, 2);
    springEnergyGradient = k*(pos-other_pos)*(length-rest_length)/length;
}

bool Connector::checkConnectivity()
{
    removed = length>=break_ratio*rest_length;
    return removed;
}

void DiscreteElasticRodsConnected::addConnector(const std::shared_ptr<Connector>& connector)
{
    connectors.push_back(connector);
}

void DiscreteElasticRodsConnected::updatePosition()
{
    Eigen::MatrixX3d prev_d3 = unitTangents(x);
    updateCenterlinePosition(); // xi+1 = xi + dt * vi
    Eigen::MatrixX3d d3 = unitTangents(x);
    updateMaterialFrame(prev_d3, d3);

    updateEdge();
    updateLength();
    updateCurvatureBinormal(d3);
}

void DiscreteElasticRodsConnected::updateConnectors()
{
    std::vector<std::shared_ptr<Connector> > new_connectors;
    for (auto& c : connectors) {
        if (!c->removed) {
            new_connectors.push_back(c);
        }
    }
    connectors = new_connectors;
}

void DiscreteElasticRodsConnected::updateVelocity()
{
    Eigen::MatrixX3d d3 = unitTangents(x);
    Eigen::VectorXd twist = getTwist(d2_ref, d3);

    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> hessian;
    if (verbose) std::cout << "- Computing Energy" << std::endl;
    computeGradientAndHessian(gradient, hessian, d3, twist);
    updateCenterlineVelocity(gradient); // vi+1 = vi + dt * F/m
    kinetic_energy = computeKineticEnergy();
    updateFrameTheta(gradient);
}

void DiscreteElasticRodsConnected::computeGradientAndHessian(Eigen::VectorXd& gradient,
        Eigen::SparseMatrix<double>& hessian,
        Eigen::MatrixX3d& d3,
        Eigen::VectorXd& twist)
{
    Eigen::VectorXd stretching_force, bending_force, twisting_force;
    stretching_force.resize(3*nv);
    stretching_force.setZero();
    bending_force.resize(3*nv);
    bending_force.setZero();
    twisting_force.resize(3*nv);
    twisting_force.setZero();

    std::tie(gradient, hessian) = createZeroGradientAndHessian();
    std::vector<Eigen::Triplet<double>>
            hessian_triplets;

    if (params.stretching_energy_enabled) {
        if (verbose) std::cout << "    - Computing Stretching Energy" << std::endl;
        stretching_energy = applyStretchingForce(gradient, hessian_triplets, d3, stretching_force);
    }

    if (params.bending_energy_enabled) {
        if (verbose) std::cout << "    - Computing Bending Energy" << std::endl;
        bending_energy = applyBendingForce(gradient, hessian_triplets, d3, bending_force);
    }

    if (params.twisting_energy_enabled) {
        if (verbose) std::cout << "    - Computing Twisting Energy" << std::endl;
        twisting_energy = applyTwistingForce(gradient, hessian_triplets, twist, twisting_force);
    }

    if (params.spring_energy_enabled) {
        if (verbose) std::cout << "    - Computing Spring Energy" << std::endl;
        twisting_energy = applySpringForce(gradient, hessian_triplets);
    }

    if (params.gravity_enabled) {
        if (verbose) std::cout << "    - Adding Gravity" << std::endl;
        potential_energy = applyGravity(gradient);
    }

    if (applyExternalForce != nullptr) {
        if (verbose) std::cout << "    - Adding External Force" << std::endl;
        applyExternalForce(gradient);
    }

    hessian.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());

    buildForceVisualization(stretching_force, bending_force, twisting_force);
}

double DiscreteElasticRodsConnected::applySpringForce(Eigen::VectorXd& gradient,
        std::vector<Eigen::Triplet<double> >& hessian_triplets)
{
    double E_s = 0.;
    for (const auto& i : connectors) {
        if (i->rod_index == index)
            gradient.segment<3>(i->node_index) += i->springEnergyGradient;
        else if (i->other_rod_index == index)
            gradient.segment<3>(i->other_node_index) += i->springEnergyGradient;
        E_s += i->springEnergy;
    }
    return E_s;
}
