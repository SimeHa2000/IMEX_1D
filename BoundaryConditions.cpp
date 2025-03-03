#include "BoundaryConditions.h"

void transmissiveBC(stateVec& consVec)
{
   consVec[0] = consVec[1];
   consVec[N+1] = consVec[N];
}