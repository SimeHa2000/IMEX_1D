#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

#include "parms.h"
#include "EoS.h"

template <typename T>
void transmissiveBC(std::vector<T>& consVec)
{
   consVec[1] = consVec[2];
   consVec[0] = consVec[1];

   consVec[N+1] = consVec[N];
   consVec[N+1] = consVec[N+1];
}

#endif

