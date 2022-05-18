#ifndef __DISCRETE_ELASTIC_RODS_VISUALIZATION_H__
#define __DISCRETE_ELASTIC_RODS_VISUALIZATION_H__
// std
#include <string>
// polyscope
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
// project
#include "DiscreteElasticRods.h"

class DiscreteElasticRodsVisualization {
public:
    polyscope::CurveNetwork* mesh = nullptr;
    std::vector<Eigen::Vector3d> vis_nodes;
    std::vector<std::array<size_t, 2>> vis_edges;

    DiscreteElasticRodsVisualization() = default;

    void initVisualization(DiscreteElasticRods& discrete_elastic_rods, std::string name)
    {
        for (int i = 0; i<discrete_elastic_rods.nv; i++) {
            vis_nodes.emplace_back(discrete_elastic_rods.x(3*i), discrete_elastic_rods.x(3*i+1),
                    discrete_elastic_rods.x(3*i+2));
        }
        for (size_t i = 1; i<vis_nodes.size(); i++) {
            std::array<size_t, 2> edge{i-1, i};
            vis_edges.push_back(edge);
        }

        mesh = polyscope::registerCurveNetwork(name, vis_nodes, vis_edges);
        mesh->setRadius(0.003);
        //mesh->setEnabled(false);
    }

    void updateVisualization(DiscreteElasticRods& discrete_elastic_rods)
    {
        for (int i = 0; i<discrete_elastic_rods.nv; i++) {
            Eigen::Vector3d node(discrete_elastic_rods.x(3*i), discrete_elastic_rods.x(3*i+1),
                    discrete_elastic_rods.x(3*i+2));
            vis_nodes[i] = node;
        }
        // update curve network
        //mesh->setEnabled(false);
        mesh->updateNodePositions(vis_nodes);
        mesh->removeAllQuantities();
        mesh->addNodeVectorQuantity("Twisting Force", discrete_elastic_rods.vis_twisting_force);
        mesh->addNodeVectorQuantity("Bending Force", discrete_elastic_rods.vis_bending_force);
        mesh->addNodeVectorQuantity("Stretching Force", discrete_elastic_rods.vis_stretching_force);
        mesh->addNodeVectorQuantity("Gradient", discrete_elastic_rods.vis_gradient);
        mesh->addEdgeVectorQuantity("d1", discrete_elastic_rods.d1);
        mesh->addEdgeVectorQuantity("d2", discrete_elastic_rods.d2);
        mesh->getQuantity("d1")->setEnabled(false);
        mesh->getQuantity("d2")->setEnabled(false);
//        Eigen::MatrixX3d node_kb;
//        node_kb.resize(discrete_elastic_rods.nv, 3);
//        node_kb.setZero();
//        for (int i = 1; i<discrete_elastic_rods.nv-1; i++)
//            node_kb.row(i) = discrete_elastic_rods.kb.row(i-1);
//        mesh->addNodeVectorQuantity("kb", node_kb);
    }
};

#endif
