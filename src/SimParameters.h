#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H
#include <cstddef>

struct SimParameters {
  double time_step;

  int newton_max_iters;
  double newton_tolerance;

  bool gravity_enabled;
  double gravity_G;

  double stretching_modulus;
  double bending_modulus;
  double twisting_modulus;
  double segment_radius;

  bool stretching_energy_enabled;
  bool bending_energy_enabled;
  bool twisting_energy_enabled;

  SimParameters()
  {
      time_step = 0.005;

      newton_max_iters = 20;
      newton_tolerance = 1e-8;

      gravity_enabled = true;
      gravity_G = -9.8;

      stretching_modulus = 1000;
      bending_modulus = 1000000;
      twisting_modulus = 10000;
      segment_radius = 0.1;

      stretching_energy_enabled = true;
      bending_energy_enabled = true;
      twisting_energy_enabled = true;
  }
};

#endif
