#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include "parms.h"

void setRiemannProblem(stateVec& cons, stateVec& prim, int argc, char** argv);

#endif // INITIALCONDITIONS_H