#ifndef __DISCRETE_ELASTIC_RODS_DRIVER_H__
#define __DISCRETE_ELASTIC_RODS_DRIVER_H__

// polyscope
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
// project
#include "DiscreteElasticRods.h"
#include "SimParameters.h"

class DiscreteElasticRodsDriver {
public:
    DiscreteElasticRods discrete_elastic_rods;
    int test = 0;

    // visualization
    polyscope::CurveNetwork* mesh = nullptr;
    std::vector<Eigen::Vector3d> vis_nodes;
    std::vector<std::array<size_t, 2>> vis_edges;

    // status
    bool running = false;

    DiscreteElasticRodsDriver() = default;

    void initVisualization(int nv, Eigen::VectorXd x_)
    {
        for (int i = 0; i<nv; i++) {
            Eigen::Vector3d node(x_(3*i), x_(3*i+1), x_(3*i+2));
            vis_nodes.push_back(node);
        }
        for (size_t i = 1; i<vis_nodes.size(); i++) {
            std::array<size_t, 2> edge{i-1, i};
            vis_edges.push_back(edge);
        }
    }

    void updateVisualization(int nv, Eigen::VectorXd x_)
    {
        for (int i = 0; i<nv; i++) {
            Eigen::Vector3d node(x_(3*i), x_(3*i+1), x_(3*i+2));
            vis_nodes[i] = node;
        }
        // update curve network
        mesh->updateNodePositions(vis_nodes);
        mesh->removeAllQuantities();
        mesh->addNodeVectorQuantity("Twisting Force", discrete_elastic_rods.vis_twisting_force);
        mesh->addNodeVectorQuantity("Bending Force", discrete_elastic_rods.vis_bending_force);
        mesh->addNodeVectorQuantity("Stretching Force", discrete_elastic_rods.vis_stretching_force);
        mesh->addNodeVectorQuantity("Gradient", discrete_elastic_rods.vis_gradient);
        mesh->addEdgeVectorQuantity("d1", discrete_elastic_rods.d1);
        mesh->addEdgeVectorQuantity("d2", discrete_elastic_rods.d2);
        mesh->getQuantity("d1")->setEnabled(true);
        mesh->getQuantity("d2")->setEnabled(true);
//        Eigen::MatrixX3d node_kb;
//        node_kb.resize(discrete_elastic_rods.nv, 3);
//        node_kb.setZero();
//        for (int i = 1; i<discrete_elastic_rods.nv-1; i++)
//            node_kb.row(i) = discrete_elastic_rods.kb.row(i-1);
//        mesh->addNodeVectorQuantity("kb", node_kb);
    }

    void initialize()
    {
        int nv = 0;
        Eigen::VectorXd x_;
        std::vector<bool> is_fixed;
        SimParameters params;
        Eigen::VectorXd theta;

        switch (test) {
        case 0:std::cout << "case 0: testing stretching and bending energy" << std::endl;
            nv = 10;
            x_.resize(nv*3);
            x_.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x_(3*i) = (-(nv >> 1)+i);
                is_fixed.emplace_back(false);
            }

            for (int i = 0; i<2; i++) {
                is_fixed[i] = true;
            }

            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = true;
            params.twisting_energy_enabled = false;
            params.gravity_enabled = true;
            break;
        case 1:std::cout << "case 1: testing twisting energy" << std::endl;
            nv = 10;
            x_.resize(nv*3);
            x_.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x_(3*i) = (-(nv >> 1)+i);
                is_fixed.emplace_back(false);
            }

            is_fixed[0] = true;
            is_fixed[1] = true;
            is_fixed[nv-2] = true;
            is_fixed[nv-1] = true;
            theta(0) = M_PI/2;
            theta(nv-2) = -M_PI/2;

            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = false;
            params.twisting_energy_enabled = true;
            params.gravity_enabled = true;
            break;
        default:std::cout << "invalid test case" << std::endl;
        }
        discrete_elastic_rods.initSimulation(nv, x_, theta, is_fixed, params);
        initVisualization(nv, x_);
        // init curve network
        mesh = polyscope::registerCurveNetwork("Discrete Elastic Rods", vis_nodes, vis_edges);
        updateVisualization(discrete_elastic_rods.nv, discrete_elastic_rods.x);
    }

    void simulateOneStep()
    {
        const unsigned int ninnersteps = 20;
        for (int inner = 0; inner<ninnersteps; inner++) {
            //std::cout << "======== step " << inner << " ========" << std::endl;
            discrete_elastic_rods.simulateOneStep();
        }
        updateVisualization(discrete_elastic_rods.nv, discrete_elastic_rods.x);
    }
};

#endif
