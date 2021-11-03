#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"

#include "DiscreteElasticRods.h"
#include "SimParameters.h"
#include "InputParser.h"

std::tuple<int, Eigen::VectorXd, Eigen::VectorXd, std::vector<bool>, SimParameters>
setTestCases(int case_number) {
    int nv = 0;
    Eigen::VectorXd x_;
    std::vector<bool> is_fixed;
    SimParameters params;
    Eigen::VectorXd theta;

    switch(case_number) {
        case 0:
            std::cout << "case 0: testing stretching energy" << std::endl;
            nv = 3;
            x_.resize(nv * 3);
            theta.resize(nv-1);
            theta.setZero();
            x_ << -1, 0, 0,
                    0, 0, 0,
                    1, 0, 0;
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(false);
            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = false;
            params.twisting_energy_enabled = false;
            params.gravity_enabled = true;
            break;
        case 1:
            std::cout << "case 0: testing bending energy" << std::endl;
            nv = 3;
            x_.resize(nv * 3);
            theta.resize(nv-1);
            theta.setZero();
            x_ << -1, 0, 0,
                    0, 0, 0,
                    1, 0, 0;
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(false);
            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = true;
            params.twisting_energy_enabled = false;

            params.gravity_enabled = true;
            break;
        case 2:
            std::cout << "case 0: testing bending energy (rge)" << std::endl;
            nv = 5;
            x_.resize(nv * 3);
            theta.resize(nv-1);
            theta.setZero();
            x_ <<   -2, 0, 0,
                    -1, 0, 0,
                    0, 0, 0,
                    1, 0, 0,
                    2, 0, 0;
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(false);
            is_fixed.emplace_back(false);
            is_fixed.emplace_back(false);
            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = true;
            params.twisting_energy_enabled = false;

            params.gravity_enabled = true;
            break;
        case 3:
            std::cout << "case 3: testing gravity" << std::endl;
            nv = 3;
            x_.resize(nv * 3);
            theta.resize(nv-1);
            theta.setZero();
            x_ << -1, 0, 0,
                    0, 0, 0,
                    1, 0, 0;
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(false);
            params.stretching_energy_enabled = false;
            params.bending_energy_enabled = false;
            params.twisting_energy_enabled = false;
            params.gravity_enabled = true;
            break;
        case 4:
            std::cout << "case 4: testing twisting" << std::endl;
            nv = 3;
            x_.resize(nv * 3);
            theta.resize(nv-1);
            theta.setZero();
            x_ << -1, 0, 0,
                    0, 0, 0,
                    1, 0, 0;
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(true);
            is_fixed.emplace_back(false);
            theta(0) = M_PI/2; //fixed

            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = false;
            params.twisting_energy_enabled = true;
            params.gravity_enabled = false;
            break;
        default:
            std::cout << "invalid case" << std::endl;
    }

    return std::make_tuple(nv, x_, theta, is_fixed, params);
}

void initVisualization(int nv, Eigen::VectorXd x_, std::vector<Eigen::Vector3d>& vis_nodes)
{
    for (int i = 0; i < nv; i++) {
        Eigen::Vector3d node(x_(3*i),x_(3*i+1),x_(3*i+2));
        vis_nodes.push_back(node);
    }
}

void buildVisualization(int nv, Eigen::VectorXd x_, std::vector<Eigen::Vector3d>& vis_nodes)
{
    for (int i = 0; i < nv; i ++) {
        Eigen::Vector3d node(x_(3*i),x_(3*i+1),x_(3*i+2));
        vis_nodes[i] = node;
    }
}

int main(int argc, char *argv[])
{
    int case_number = 0;

    InputParser input(argc, argv);
    if (input.cmdOptionExists("-test")) {
        std::string input_test_number = input.getCmdOption("-test");
        if (!input_test_number.empty()) case_number = std::stoi(input_test_number);
    }

    int nv = 0;
    Eigen::VectorXd x_;
    Eigen::VectorXd theta;
    std::vector<bool> is_fixed;
    SimParameters params;
    std::tie(nv, x_, theta, is_fixed, params) = setTestCases(case_number);

    polyscope::init();

    std::vector<Eigen::Vector3d> vis_nodes;
    std::vector<std::array<size_t, 2>> vis_edges;

    // convert x to 3D nodes
    initVisualization(nv, x_, vis_nodes);
    for (size_t i = 1; i < vis_nodes.size(); i++) {
        std::array<size_t, 2> edge{i-1, i};
        vis_edges.push_back(edge);
    }

    // Add the curve network
    polyscope::CurveNetwork* mesh = polyscope::registerCurveNetwork("Mesh 0", vis_nodes, vis_edges);

    // solver iteration steps setting
    int noutersteps = 10;    // number of minimal surface optimization steps
    int ninnersteps = 20;  // number of newton's method steps

    DiscreteElasticRods discrete_elastic_rods;
    discrete_elastic_rods.initSimulation(nv, x_, theta, is_fixed, params);

    for (int outer = 0; outer < noutersteps; outer++)
    {
        for (int inner = 0; inner < ninnersteps; inner++)
        {
            std::cout << "======== step " << inner << " ========" << std::endl;
            discrete_elastic_rods.simulateOneStep();
        }
        buildVisualization(nv, discrete_elastic_rods.x, vis_nodes);

        // add the mesh at this optimization step tp viewer
        mesh = polyscope::registerCurveNetwork("Mesh " + std::to_string(outer+1), vis_nodes, vis_edges);

        mesh->addNodeVectorQuantity("Twisting Force", discrete_elastic_rods.vis_twisting_force);
        mesh->addNodeVectorQuantity("Bending Force", discrete_elastic_rods.vis_bending_force);
        mesh->addNodeVectorQuantity("Stretching Force", discrete_elastic_rods.vis_stretching_force);
        mesh->addNodeVectorQuantity("Gradient", discrete_elastic_rods.vis_gradient);
        mesh->addEdgeVectorQuantity("d1", discrete_elastic_rods.d1);
        mesh->addEdgeVectorQuantity("d2", discrete_elastic_rods.d2);
        Eigen::MatrixX3d node_kb;
        node_kb.resize(discrete_elastic_rods.nv, 3);
        node_kb.setZero();
        for (int i = 1; i < discrete_elastic_rods.nv-1; i++)
            node_kb.row(i) = discrete_elastic_rods.kb.row(i-1);
        mesh->addNodeVectorQuantity("kb", node_kb);
    }

    // visualize
    polyscope::show();
}