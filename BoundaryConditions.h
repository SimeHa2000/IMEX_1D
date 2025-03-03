#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

#include "parms.h"

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

void transmissiveBC(stateVec& consVec);

#endif

