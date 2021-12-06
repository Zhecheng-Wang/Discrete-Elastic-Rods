#ifndef __CONNECTOR_H__
#define __CONNECTOR_H__

// std
#include <vector>
#include <tuple>
#include <functional>
// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
// project
#include "DiscreteElasticRods.h"

class DiscreteElasticRodsConnected;
typedef DiscreteElasticRodsConnected DERC;

class Connector {
public:
    int index;
    int removed = false;

    int rod_index;
    int node_index;
    double segment_ratio;

    int other_rod_index;
    int other_node_index;
    double other_segment_ratio;

    double k;
    double length;
    double rest_length;
    double break_ratio = 2.;
    double springEnergy;
    Eigen::Vector3d springEnergyGradient;

    Connector() = default;

    void initConnector(std::vector<std::shared_ptr<DERC> > rods,
            int _rod_index,
            std::pair<int, double>& node_index_segment_ratio_pair,
            int _other_rod_index,
            std::pair<int, double>& other_node_index_segment_ratio_pair,
            double _k);

    void computeConnectorEnergy(std::vector<std::shared_ptr<DERC> > rods);

    bool checkConnectivity();
};

class DiscreteElasticRodsConnected : public DiscreteElasticRods {
public:
    int index;
    std::vector<std::shared_ptr<Connector> > connectors;

    DiscreteElasticRodsConnected() = default;

    void addConnector(const std::shared_ptr<Connector>& connector);

    void computeGradientAndHessian(Eigen::VectorXd& gradient,
            Eigen::SparseMatrix<double>& hessian,
            Eigen::MatrixX3d& d3,
            Eigen::VectorXd& twist) override;

    void updatePosition();

    void updateConnectors();

    void updateVelocity();

    double applySpringForce(Eigen::VectorXd& gradient,
            std::vector<Eigen::Triplet<double> >& hessian_triplets);
};

#endif
