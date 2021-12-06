// std
#include <iostream>
// polyscope
#include "polyscope/polyscope.h"
// project
#include "InputParser.h"
#include "DiscreteElasticRodsDriver.h"
//#include "DiscreteElasticRodsAMDriver.h"

DiscreteElasticRodsDriver discrete_elastic_rods_driver;
//DiscreteElasticRodsAMDriver discrete_elastic_rods_am_driver;

void discreteElasticRodsDriverCallback()
{
    ImGui::Text("Status: %s", (discrete_elastic_rods_driver.running) ? "Running" : "Stopped");
    ImGui::Text("Stretching Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.stretching_energy);
    ImGui::Text("Bending Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.bending_energy);
    ImGui::Text("Twisting Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.twisting_energy);
    ImGui::Text("Potential Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.potential_energy);
    ImGui::Text("Kinetic Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.kinetic_energy);
    ImGui::Text("Total Energy = %f", discrete_elastic_rods_driver.discrete_elastic_rods.stretching_energy +
                                        discrete_elastic_rods_driver.discrete_elastic_rods.bending_energy +
                                        discrete_elastic_rods_driver.discrete_elastic_rods.twisting_energy +
                                        discrete_elastic_rods_driver.discrete_elastic_rods.potential_energy +
                                        discrete_elastic_rods_driver.discrete_elastic_rods.kinetic_energy);
    ImGui::Checkbox("stretching", &discrete_elastic_rods_driver.discrete_elastic_rods.params.stretching_energy_enabled);
    ImGui::SameLine();
    ImGui::Checkbox("bending", &discrete_elastic_rods_driver.discrete_elastic_rods.params.bending_energy_enabled);
    ImGui::SameLine();
    ImGui::Checkbox("twisting", &discrete_elastic_rods_driver.discrete_elastic_rods.params.twisting_energy_enabled);
    ImGui::SameLine();
    ImGui::Checkbox("gravity", &discrete_elastic_rods_driver.discrete_elastic_rods.params.gravity_enabled);
    ImGui::InputDouble("time step dt", &discrete_elastic_rods_driver.discrete_elastic_rods.params.time_step);
    if (ImGui::Button("Run/Stop Simulation"))
        discrete_elastic_rods_driver.running = !discrete_elastic_rods_driver.running;
    ImGui::SameLine();
    if (ImGui::Button("Run One Time Step Simulation")) discrete_elastic_rods_driver.simulateOneStep();
    ImGui::SameLine();
    if (ImGui::Button("Reset Simulation")) discrete_elastic_rods_driver.resetSimulation();
    if (discrete_elastic_rods_driver.running) discrete_elastic_rods_driver.simulateOneStep();
}

//void discreteElasticRodsAMDriverCallback()
//{
//    ImGui::Text("Status: %s", (discrete_elastic_rods_am_driver.running) ? "Running" : "Stopped");
//    if (ImGui::Button("Run/Stop Simulation"))
//        discrete_elastic_rods_am_driver.running = !discrete_elastic_rods_am_driver.running;
//    ImGui::SameLine();
//    if (ImGui::Button("Run One Time Step Simulation")) discrete_elastic_rods_am_driver.simulateOneStep();
//    if (discrete_elastic_rods_am_driver.running) discrete_elastic_rods_am_driver.simulateOneStep();
//}

int main(int argc, char* argv[])
{
    InputParser input(argc, argv);
    int driver = input.getIntegerCmdOption("-driver");
    int test = input.getIntegerCmdOption("-test");

    // polyscope setup
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
    switch (driver) {
    case 0: {
        polyscope::options::programName = "Discrete Elastic Rods Simulation";
        polyscope::options::printPrefix = "[DER] ";
    }
        break;
    case 1: {
        polyscope::options::programName = "Discrete Elastic Rods Based Additive Manufactured Parts Simulation";
        polyscope::options::printPrefix = "[DER-AM] ";
    }
        break;
    default: {
        std::cerr << "[ERROR] Please Select A Valid Driver.";
    }
    }

    // driver setup
    polyscope::init();
    switch (driver) {
    case 0: {
        discrete_elastic_rods_driver.test = test;
        discrete_elastic_rods_driver.initialize();
        polyscope::state::userCallback = discreteElasticRodsDriverCallback;
    }
        break;
    case 1: {
//        discrete_elastic_rods_am_driver.test = test;
//        discrete_elastic_rods_am_driver.initialize();
//        polyscope::state::userCallback = discreteElasticRodsAMDriverCallback;
    }
        break;
    default: {
        std::cerr << "[ERROR] Please Select A Valid Driver.";
    }
    }
    polyscope::show();
}
