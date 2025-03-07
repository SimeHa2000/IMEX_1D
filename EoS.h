#ifndef EOS_H
#define EOS_H

#include "parms.h"

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

extern double Gamma;
extern int N;

void primToCons(state &prim, state &cons);
void consToPrim(state &prim, state &cons);
double getPressure(state &cons);
double getKineticEnergy(state& cons);
double getEnthalpy(state& cons);

#endif
