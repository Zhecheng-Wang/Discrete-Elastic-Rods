#ifndef __DISCRETE_ELASTIC_RODS_AM_H__
#define __DISCRETE_ELASTIC_RODS_AM_H__

#include "DiscreteElasticRodsConnected.h"

class DiscreteElasticRodsAM {
public:
    std::vector<std::shared_ptr<DERC> > rods;
    std::vector<std::shared_ptr<Connector> > connectors;

    DiscreteElasticRodsAM();

    void initSimulation(std::vector<std::shared_ptr<DERC> > _rods);

    void simulateOneStep();
};

#endif
