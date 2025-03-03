#ifndef IMEX_H
#define IMEX_H

#include "parms.h"

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <Eigen/Sparse>

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

void compute_dt(stateVec& consVec, double& dt);
void explicitFluxes(stateVec& consVec, stateVec& expFluxes);
void implicitFLuxes(stateVec& consVec, stateVec& impFluxes);
void explicitUpdate(stateVec& consVec, stateVec& expFluxes, double dt);

#endif