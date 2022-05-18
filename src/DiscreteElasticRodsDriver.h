#ifndef __DISCRETE_ELASTIC_RODS_DRIVER_H__
#define __DISCRETE_ELASTIC_RODS_DRIVER_H__

// polyscope
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
// project
#include "DiscreteElasticRods.h"
#include "DiscreteElasticRodsVisualization.h"
#include "SimParameters.h"

typedef std::tuple<int, Eigen::VectorXd, Eigen::VectorXd, std::vector<bool>, SimParameters> RodInitConfiguration;

class DiscreteElasticRodsDriver {
public:
    DiscreteElasticRods discrete_elastic_rods;
    int test = 0;

    // visualization
    DiscreteElasticRodsVisualization discrete_elastic_rods_visualization;

    // status
    bool running = false;
    RodInitConfiguration init_state;

    DiscreteElasticRodsDriver() = default;

    void initialize()
    {
        int nv = 0;
        Eigen::VectorXd x;
        std::vector<bool> is_fixed;
        SimParameters params;
        Eigen::VectorXd theta;

        switch (test) {
        case 0:
            std::cout << "case 0: testing stretching and bending energy" << std::endl;
            nv = 4;
            x.resize(nv*3);
            x.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x(3*i) = (-(nv >> 1)+i);
                is_fixed.emplace_back(false);
            }

            is_fixed[0] = true;

            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = true;
            params.twisting_energy_enabled = false;
            params.gravity_enabled = true;
            break;
        case 1:
            std::cout << "case 1: testing twisting energy" << std::endl;
            nv = 10;
            x.resize(nv*3);
            x.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x(3*i) = (-(nv >> 1)+i);
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
        case 2:
            std::cout << "case 2: testing external force" << std::endl;
            nv = 3;
            x.resize(nv*3);
            x.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x(3*i) = (-(nv >> 1)+i);
                is_fixed.emplace_back(false);
            }

            is_fixed[0] = true;

            params.stretching_energy_enabled = true;
            params.bending_energy_enabled = false;
            params.twisting_energy_enabled = false;
            params.gravity_enabled = true;

            discrete_elastic_rods.applyExternalForce = [this](Eigen::VectorXd& gradient) {
                gradient(7) -= 9.0;
            };
            break;
        default:
            std::cout << "invalid test case" << std::endl;
        }
        init_state = std::make_tuple(nv, x, theta, is_fixed, params);
        discrete_elastic_rods.initSimulation(nv, x, theta, is_fixed, params);
        discrete_elastic_rods_visualization.initVisualization(discrete_elastic_rods, "Discrete Elastic Rods");
        discrete_elastic_rods_visualization.updateVisualization(discrete_elastic_rods);
    }

    void simulateOneStep()
    {
        const unsigned int ninnersteps = 20;
        for (int inner = 0; inner<ninnersteps; inner++) {
            //std::cout << "======== step " << inner << " ========" << std::endl;
            discrete_elastic_rods.simulateOneStep();
        }
        discrete_elastic_rods_visualization.updateVisualization(discrete_elastic_rods);
    }

    void resetSimulation()
    {
        discrete_elastic_rods.initSimulation(std::get<0>(init_state), std::get<1>(init_state), std::get<2>(init_state),
                std::get<3>(init_state), std::get<4>(init_state));
        discrete_elastic_rods_visualization.updateVisualization(discrete_elastic_rods);
    }
};

#endif
