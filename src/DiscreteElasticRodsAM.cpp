#include "DiscreteElasticRodsAM.h"

DiscreteElasticRodsAM::DiscreteElasticRodsAM() = default;

void DiscreteElasticRodsAM::initSimulation(std::vector<std::shared_ptr<DERC> > _rods)
{
    rods = _rods;
    int nrods = rods.size();
    for (int i = 0; i < nrods; i++) {
        DERC& rod = *rods[i];
        int nv = rod.nv;
        for (int j = i; j < nrods; j++) {
            DERC& other_rod = *rods[j];
            int other_nv = other_rod.nv;
            for (int j )
        }
    }
}

void DiscreteElasticRodsAM::simulateOneStep()
{
    for (const auto& der : rods) {
        der->updatePosition();
    }
    for (const auto& c : connectors) {
        c->computeConnectorEnergy(rods);
    }
    for (const auto& der: rods) {
        der->updateVelocity();
    }
}
