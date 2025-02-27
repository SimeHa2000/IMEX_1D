#include "IMEX.h"
#include "parms.h"
#include "EoS.h" 

extern double x0;
extern double x1;
extern double epsilon; 
extern double Gamma;
extern double CFL;
extern int N;
double dx = (x1 - x0) / N;

void compute_dt(stateVec& consVec, double& dt){
    
    double wave_speed;

    for(int i = 0; i < N; i++){
        double rho = consVec[i][0];
        double velocity = consVec[i][1]/rho;
        double pressure = getPressure(consVec[i]);

        double max_speed = 2.0 * abs(velocity);

        if(max_speed > wave_speed){
            wave_speed = max_speed;
        }
    }

    dt = CFL * dx / wave_speed;
}


void explicitFluxes(stateVec& consVec, stateVec& expFluxes){
    // explicit flux from the Euler equations using VT splitting

    for(int i=0; i < consVec.size(); i++){
        state cons = consVec[i];
        state expFlux;

        double rho = cons[0];
        double mom = cons[1];
        double E = cons[2];

        double k = getKineticEnergy(cons);  // get kinetic energy 

        expFlux[0] = mom; 
        expFlux[1] = (mom * mom)/rho;
        expFlux[2] = k * mom/rho;

        expFluxes[i] = expFlux;
    }
}

void implicitFLuxes(stateVec& consVec, stateVec& impFluxes){
    // implicit flux from the Euler equations using VT splitting
    for(int i=0; i < consVec.size(); i++){
        state cons = consVec[i];
        state impFlux;
    
        double rho = cons[0];
        double mom = cons[1];
        double E = cons[2];

        double k = getKineticEnergy(cons);  // get kinetic energy 
        double p = getPressure(cons); // get pressure
        double h = getEnthalpy(cons); // get enthalpy


        impFlux[0] = 0.0; 
        impFlux[1] = p/epsilon;
        impFlux[2] = h * mom/rho;

        impFluxes[i] = impFlux;
        
    }
}

void explicitUpdate(stateVec& F_exp, stateVec& consVec){

    stateVec expFluxes(N+nGhost);
    stateVec fPluaHalf(N+nGhost);
    stateVec fMinusHalf(N+nGhost);

    double aSpeed, dt;

    explicitFluxes(consVec, expFluxes);
    compute_dt(consVec, dt);

    for(int i = 1; i <= N; i++){

        aSpeed = std::max(std::fabs(consVec[i][1]/consVec[i][0]) , std::fabs(consVec[i+1][1]/consVec[i+1][0]));

        for(int var = 0; var < nVars; var++){
            fPluaHalf[i][var] = 0.5 * (expFluxes[i+1][var] + expFluxes[i][var]) - 0.5 * std::fabs(aSpeed) * (consVec[i+1][var] - consVec[i][var]);
            fMinusHalf[i][var] = 0.5 * (expFluxes[i][var] + expFluxes[i-1][var]) - 0.5 * std::fabs(aSpeed) * (consVec[i][var] - consVec[i-1][var]);
        }
    }

    //Final explicit update
    for(int i = 1; i <= N; i++){

        for(int var = 0; var < nVars; var++){
            F_exp[i][var] = consVec[i][var] - dt/dx * (fPluaHalf[i][var] - fMinusHalf[i][var]);
        }
    }
}

void implicitUpdate(stateVec& F_imp, stateVec& consVec){

}

void IMEXupdate(
    
)


