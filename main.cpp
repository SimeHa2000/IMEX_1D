#include "parms.h"
#include "EoS.h"
#include "IMEX.h"
#include "BoundaryConditions.h"

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

extern double Gamma;
extern int N;
extern double epsilon;
extern double CFL;
extern double x0;
extern double x1;
extern int nGhost;

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;


void setRiemannProblem(stateVec& cons, stateVec& prim, int argc, char** argv){
    
    state stateL;
    state stateR;

    readStates(stateL, stateR, argc, argv);

    // asser that N is even
    assert(N % 2 == 0);

    for(int i = 1; i <= N; i++)
    {
        if(i < N/2){
            cons[i] = stateL;
        }
        
        else
        cons[i] = stateR;
        
    consToPrim(prim[i], cons[i], epsilon);
    }
}

int main(int argc, char** argv)
{
    // set up variables

    readParams(argc, argv);
    int nGhost = 2;

    stateVec consVec(N + nGhost);
    stateVec primVec(N + nGhost);

    setRiemannProblem(consVec, primVec, argc, argv);
    transmissiveBC(consVec);
}
