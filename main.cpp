#include "parms.h"
#include "EoS.h"
#include "IMEX.h"

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


void variableSetup(state cons, state prim)
{
    cons[0] = 0.0;
    cons[1] = 0.0;
    cons[2] = 0.0;

    prim[0] = 0.0;
    prim[1] = 0.0;
    prim[2] = 0.0;

}

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

void transmissiveBC(stateVec& consVec)
{
   consVec[0] = consVec[1];
   consVec[N+1] = consVec[N];

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


    // while (t <= t_stop)
    // {
    //     compute_dt(cons, prim);

    //     state consImp;
    //     state consExp;

    //     Rusanov_adv(consExp);

    //     Picard(consImp);
    // }
}
