#ifndef __DISCRETE_ELASTIC_RODS_DRIVER_AM_H__
#define __DISCRETE_ELASTIC_RODS_DRIVER_AM_H__

// project
#include "DiscreteElasticRods.h"
#include "DiscreteElasticRodsAM.h"
#include "DiscreteElasticRodsVisualization.h"
#include "SimParameters.h"

typedef std::tuple<int, Eigen::VectorXd, Eigen::VectorXd, std::vector<bool>, SimParameters> RodInitConfiguration;

class DiscreteElasticRodsAMDriver {
public:
    DiscreteElasticRodsAM discrete_elastic_rods_am;
    int test = 0;

    // visualization
    std::vector<DiscreteElasticRodsVisualization> rods_visualization;

    // status
    bool running = false;
    std::vector<RodInitConfiguration> init_states;

    DiscreteElasticRodsDriver() = default;

    void initialize()
    {
        std::vector<int> nv;
        std::vector<Eigen::VectorXd> x;
        std::vector<bool> is_fixed;
        std::vector<SimParameters> params;
        std::vector<Eigen::VectorXd> theta;

        switch (test) {
        case 0:std::cout << "case 0: testing two rods" << std::endl;
            nv = 10;
            x.resize(nv*3);
            x.setZero();
            theta.resize(nv-1);
            theta.setZero();

            for (int i = 0; i<nv; i++) {
                x(3*i) = (-(nv >> 1)+i);
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
        default:std::cout << "invalid test case" << std::endl;
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
        int nrods = discrete_elastic_rods_am.rods.size();
        for (int i = 0; i < nrods; i++) {
            DiscreteElasticRods& discrete_elastic_rods = discrete_elastic_rods_am.rods[i];
            RodInitConfiguration init_state = init_states[i];
            discrete_elastic_rods.initSimulation(std::get<0>(init_state), std::get<1>(init_state), std::get<2>(init_state), std::get<3>(init_state), std::get<4>(init_state));
            rods_visualization[i].updateVisualization(discrete_elastic_rods);
        }
    }
};

#endif
