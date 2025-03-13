#include "BoundaryConditions.h"
#include "EoS.h"
#include "IMEX.h"
#include "InitialConditions.h"
#include "parms.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

extern double Gamma;
extern int    N;
extern double epsilon;
extern double CFL;
extern double x0;
extern double x1;
extern int    nGhost;

int main(int argc, char **argv)
{
    // set up variables

    readParams(argc, argv);
    int    nGhost = 2;
    double t      = 0.0;
    double dt;
    double dx   = (x1 - x0) / N;
    int    iter = 0;

    stateVec consVec(N + 2*nGhost);
    stateVec primVec(N + 2*nGhost);

    setRiemannProblem(consVec, primVec, argc, argv);
    transmissiveBC(consVec, nGhost);

    while (t < t_stop)
    {
        compute_dt(consVec, dt, iter);

        std::cout << "iteration " << iter + 1 << ",  time = " << t << ",  dt = " << dt
                  << std::endl;

        IMEXupdate(consVec, dt);

        t += dt;

        iter++;

        std::string   fileName = "output/ToroTest1/ToroTest1_100.csv";
        std::ofstream outFile(fileName);

        for (int i = nGhost; i < N + nGhost; i++)
        {
            double x = x0 + (i - nGhost+0.5) * dx;
            consToPrim(primVec[i], consVec[i]);
            outFile << x << " , " << primVec[i][0] << " , " << primVec[i][1] << " , "
                    << primVec[i][2] << std::endl;

            //non-conservative update 

        }
        outFile.close();
    }
    return 0;
}
