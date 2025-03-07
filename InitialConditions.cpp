
#include "InitialConditions.h"
#include "EoS.h"
#include "parms.h"

void setRiemannProblem(stateVec &cons, stateVec &prim, int argc, char **argv)
{

    state stateL;
    state stateR;

    readStates(stateL, stateR, argc, argv);

    // asser that N is even
    assert(N % 2 == 0);

    for (int i = 1; i <= N; i++)
    {
        if (i < N / 2)
        {
            prim[i] = stateL;
        }

        else
            prim[i] = stateR;

        primToCons(prim[i], cons[i]);
    }
}
