#include "EoS.h"
#include "parms.h"

extern double Gamma;
extern int N;
extern double epsilon;
extern double CFL;
extern double x0;
extern double x1;
extern int nGhost;

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

void primToCons(state &prim, state &cons, double &epsilon)
{
    double rho      = prim[0];
    double velocity = prim[1];
    double pressure = prim[2];

    cons[0] = rho;
    cons[1] = rho * velocity;
    cons[2]
        = (pressure / (Gamma - 1)) + 0.5 * epsilon * velocity * velocity * rho;
}

void consToPrim(state &prim, state &cons, double &epsilon)
{
    double rho       = cons[0];
    double mom       = cons[1];
    double totEnergy = cons[2];

    prim[0] = rho;
    prim[1] = mom / rho;
    prim[2] = (Gamma - 1) * (totEnergy - (0.5 * epsilon * (mom * mom) / rho));
}

double getPressure(state &cons)
{
    double rho       = cons[0];
    double mom       = cons[1];
    double totEnergy = cons[2];

    double pressure = (Gamma - 1) * (totEnergy - (0.5 * epsilon * (mom * mom) / rho));
    return pressure;
}

double getKineticEnergy(state& cons){

    double rho = cons[0];
    double mom = cons[1];
    double E = cons[2];

    double k = epsilon * 0.5 * mom * mom / rho;

    return k;
}

double getEnthalpy(state& cons){

    double rho = cons[0];
    double mom = cons[1];
    double E = cons[2];

    double p = getPressure(cons);
    double h = p / Gamma;

    return h;
}