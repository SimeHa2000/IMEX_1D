#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "EoS.h"
#include "parms.h"

template <typename T> void transmissiveBC(std::vector<T> &consVec, int nGhost)
{
    if (nGhost == 1)
    {
        consVec[0] = consVec[1];

        consVec[N + nGhost] = consVec[N + nGhost - 1];
    }

    else if (nGhost == 2)
    {
        consVec[1] = consVec[2];
        consVec[0] = consVec[1];

        consVec[N + nGhost]     = consVec[N + nGhost - 1];
        consVec[N + nGhost + 1] = consVec[N + nGhost];
    }

    else
    {
        std::cout << "this amount of ghost cells is not supported. Enter 1 or 2."
                  << std::endl;
    }
}

#endif
