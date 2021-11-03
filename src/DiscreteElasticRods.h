#ifndef __DISCRETE_ELASTIC_RODS_H__
#define __DISCRETE_ELASTIC_RODS_H__

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>
#include "SimParameters.h"

class DiscreteElasticRods
{
    public:
    Eigen::VectorXd x_ref;
    Eigen::VectorXd x;
    std::vector<bool> is_fixed;
    Eigen::VectorXd v;
    Eigen::VectorXd omega;
    Eigen::VectorXd e;
    Eigen::VectorXd length_rest;
    Eigen::VectorXd l_ref;
    Eigen::VectorXd length;
    Eigen::MatrixX3d kb;
    Eigen::MatrixX2d kappa_ref;
    Eigen::VectorXd twist_rest;
    int nv;

    // time-parallel frame
    Eigen::MatrixX3d d1_ref;
    Eigen::MatrixX3d d2_ref;
    Eigen::MatrixX3d d1;
    Eigen::MatrixX3d d2;
    Eigen::VectorXd theta;

    // simulation parameters
    SimParameters params;
    bool verbose = true;

    // visualization
    std::vector<Eigen::Vector3d> vis_gradient;
    std::vector<Eigen::Vector3d> vis_stretching_force;
    std::vector<Eigen::Vector3d> vis_bending_force;
    std::vector<Eigen::Vector3d> vis_twisting_force;

    DiscreteElasticRods();

    void initSimulation(int nv_, Eigen::VectorXd x_, Eigen::VectorXd theta_, std::vector<bool> is_fixed_, SimParameters params_);

    void simulateOneStep();

    void updateCenterlinePosition(void);

    void updateMaterialFrame(Eigen::MatrixX3d prev_d3, Eigen::MatrixX3d d3);

    std::tuple<Eigen::MatrixXd, Eigen::SparseMatrix<double> > createZeroGradientAndHessian();

    void computeGradientAndHessian(Eigen::VectorXd& gradient,
                                   Eigen::SparseMatrix<double>& hessian,
                                   Eigen::MatrixX3d& d3,
                                   Eigen::VectorXd& twist);

    void updateCenterlineVelocity(Eigen::VectorXd& gradient);

    void updateFrameTheta(Eigen::VectorXd& gradient);

    void updateEdge();

    void updateLength();

    void updateCurvatureBinormal(Eigen::MatrixX3d d3);

    Eigen::MatrixX3d unitTangents(const Eigen::VectorXd x_);

    Eigen::Vector3d parallelTransport(Eigen::Vector3d v, Eigen::Vector3d r1, Eigen::Vector3d r2);

    double applyStretchingForce(Eigen::VectorXd& gradient,
                                std::vector<Eigen::Triplet<double>>& hessian,
                                Eigen::MatrixX3d& d3,
                                Eigen::VectorXd& stretching_force);

    void getMaterialCurvature(Eigen::MatrixX2d& kappa);

    Eigen::VectorXd getTwist(Eigen::MatrixX3d& d2, Eigen::MatrixX3d& d3);

    void getVoronoiLength(Eigen::VectorXd& l);

    Eigen::Matrix3d getCrossMatrix(Eigen::Vector3d v);

    double applyBendingForce(Eigen::VectorXd& gradient,
                             std::vector<Eigen::Triplet<double>>& hessian,
                             Eigen::MatrixX3d& d3,
                             Eigen::VectorXd& bending_force);

    double applyTwistingForce(Eigen::VectorXd& gradient,
                              std::vector<Eigen::Triplet<double>>& hessian,
                              Eigen::VectorXd& twist,
                              Eigen::VectorXd& twisting_force);

    double applyGravity(Eigen::VectorXd& gradient);

    void buildVisualization(Eigen::VectorXd& gradient);

    void buildForceVisualization(Eigen::VectorXd& stretching_force,
                                 Eigen::VectorXd& bending_force,
                                 Eigen::VectorXd& twisting_force);

    protected:
    

    private:
};

#endif
