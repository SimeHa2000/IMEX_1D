#include "IMEX.h"
#include "EoS.h"
#include "parms.h"
#include <array>
#include <set>
#include <vector>

extern double x0;
extern double x1;
extern double epsilon;
extern double Gamma;
extern double CFL;
extern int    N;
double        dx = (x1 - x0) / N;

void compute_dt(stateVec &consVec, double &dt)
{

    double wave_speed;

    for (int i = 0; i < N; i++)
    {
        double rho      = consVec[i][0];
        double velocity = consVec[i][1] / rho;
        double pressure = getPressure(consVec[i]);

        double max_speed = 2.0 * abs(velocity);

        if (max_speed > wave_speed)
        {
            wave_speed = max_speed;
        }
    }

    dt = CFL * dx / wave_speed;
}

void explicitFluxes(stateVec &consVec, stateVec &expFluxes)
{
    // explicit flux from the Euler equations using VT splitting

    for (int i = 0; i < consVec.size(); i++)
    {
        state cons = consVec[i];
        state expFlux;

        double rho = cons[0];
        double mom = cons[1];
        double E   = cons[2];

        double k = getKineticEnergy(cons); // get kinetic energy

        expFlux[0] = mom;
        expFlux[1] = (mom * mom) / rho;
        expFlux[2] = k * mom / rho;

        expFluxes[i] = expFlux;
    }
}

void implicitFLuxes(stateVec &consVec, stateVec &impFluxes)
{
    // implicit flux from the Euler equations using VT splitting
    for (int i = 0; i < consVec.size(); i++)
    {
        state cons = consVec[i];
        state impFlux;

        double rho = cons[0];
        double mom = cons[1];
        double E   = cons[2];

        double k = getKineticEnergy(cons); // get kinetic energy
        double p = getPressure(cons);      // get pressure
        double h = getEnthalpy(cons);      // get enthalpy

        impFlux[0] = 0.0;
        impFlux[1] = p / epsilon;
        impFlux[2] = h * mom / rho;

        impFluxes[i] = impFlux;
    }
}

void explicitUpdate(stateVec &F_exp, stateVec &consVec)
{

    stateVec expFluxes(N + nGhost);
    stateVec fPluaHalf(N + nGhost);
    stateVec fMinusHalf(N + nGhost);

    double aSpeed, dt;

    explicitFluxes(consVec, expFluxes);
    compute_dt(consVec, dt);

    for (int i = 1; i <= N; i++)
    {

        aSpeed = std::max(std::fabs(consVec[i][1] / consVec[i][0]),
                          std::fabs(consVec[i + 1][1] / consVec[i + 1][0]));

        for (int var = 0; var < nVars; var++)
        {
            fPluaHalf[i][var]
                = 0.5 * (expFluxes[i + 1][var] + expFluxes[i][var])
                  - 0.5 * std::fabs(aSpeed) * (consVec[i + 1][var] - consVec[i][var]);
            fMinusHalf[i][var]
                = 0.5 * (expFluxes[i][var] + expFluxes[i - 1][var])
                  - 0.5 * std::fabs(aSpeed) * (consVec[i][var] - consVec[i - 1][var]);
        }
    }

    // Final explicit update
    for (int i = 1; i <= N; i++)
    {

        for (int var = 0; var < nVars; var++)
        {
            F_exp[i][var]
                = consVec[i][var] - dt / dx * (fPluaHalf[i][var] - fMinusHalf[i][var]);
        }
    }
}

void initialisePicardVars(std::vector<double> &enthalpies,
                          std::vector<double> &hHalf_tilde,
                          std::vector<double> &initKineticEnergy,
                          std::vector<double> &initPressure,
                          std::vector<double> &initMomentum, stateVec &consVec,
                          stateVec &explicitConsVec)
{
    for (int i = 0; i < enthalpies.size(); i++)
    {
        enthalpies[i] = getEnthalpy(explicitConsVec[i]);
    }

    for (int i = 0; i < hHalf_tilde.size(); i++)
    {
        hHalf_tilde[i]
            = 0.5
              * (((enthalpies[i] / explicitConsVec[i][0]) * explicitConsVec[i][1]
                  + (enthalpies[i + 1] / explicitConsVec[i + 1][0])
                        * explicitConsVec[i + 1][1]))
              / (explicitConsVec[i][1] + explicitConsVec[i + 1][1]);
    }

    for (int i = 0; i < initKineticEnergy.size(); i++)
    {
        initKineticEnergy[i] = epsilon / 2 * explicitConsVec[i][0]
                               * (explicitConsVec[i][1] / explicitConsVec[i][0])
                               * (explicitConsVec[i][1] / explicitConsVec[i][0]);
    }

    for (int i = 0; i < initPressure.size(); i++)
    {
        // initPressure[i] = (Gamma - 1.0) * consVec[i][2]
        //                   - kStar[i]; // TODO: checks this - kStar or k

        initPressure[i] = (Gamma - 1.0) * consVec[i][2] - getKineticEnergy(consVec[i]);
    }

    for (int i = 0; i < initMomentum.size(); i++)
    {
        initMomentum[i] = explicitConsVec[i][1]; // Is this correct?
    }
}
void setbVector(Eigen::VectorXd &b, stateVec &consVec, stateVec &explicitConsVec,
                std::vector<double> &hHalf_tilde, std::vector<double> &initKineticEnergy,
                double dt)
{

    b(0) = explicitConsVec[0][2] - initKineticEnergy[0]
           - 0.5 * dt / dx
                 * (hHalf_tilde[0] * explicitConsVec[1][1]
                    - hHalf_tilde[0] * explicitConsVec[0][1]);

    for (int i = 1; i < N - 1; i++)
    {
        b(i) = explicitConsVec[i][2] - initKineticEnergy[i]
               - 0.5 * dt / dx
                     * (hHalf_tilde[i] * explicitConsVec[i + 1][1]
                        + (hHalf_tilde[i] - hHalf_tilde[i - 1]) * explicitConsVec[i][1]
                        - hHalf_tilde[i - 1] * explicitConsVec[i - 1][1]);
    }

    b(N) = explicitConsVec[N][2] - initKineticEnergy[N]
           - 0.5 * dt / dx
                 * (hHalf_tilde[N] * explicitConsVec[N][1]
                    - hHalf_tilde[N - 1] * explicitConsVec[N - 1][1]);
}

void setTMatrix(Eigen::SparseMatrix<double> &T, std::vector<double> &hHalf_tilde,
                double dt)
{

    // setting up matrix T using eigen
    T.reserve(Eigen::VectorXi::Constant(N, 5));

    // Manually account for boundary conditions for fist and last two rows
    // Transmissive case
    // Row 1
    T.insert(0, 0) = -0.25 * dt * dt / dx * hHalf_tilde[0]
                     + 0.25 * dt * dt / dx * (hHalf_tilde[0] - hHalf_tilde[0])
                     + epsilon * dx / (Gamma - 1.0)
                     + 0.25 * dt * dt / dx * (hHalf_tilde[0] - hHalf_tilde[0]);
    T.insert(0, 1) = -0.25 * dt * dt / dx * (hHalf_tilde[0] - hHalf_tilde[0]);
    T.insert(0, 2) = -0.25 * dt * dt / dx * hHalf_tilde[0];

    // Row 2
    T.insert(1, 0) = -0.25 * dt * dt / dx * hHalf_tilde[0]
                     - 0.25 * dt * dt / dx * (hHalf_tilde[1] - hHalf_tilde[0]);
    T.insert(1, 1) = epsilon * dx / (Gamma - 1.0)
                     + 0.25 * dt * dt / dx * (hHalf_tilde[1] - hHalf_tilde[0]);
    T.insert(1, 2) = -0.25 * dt * dt / dx * (hHalf_tilde[1] - hHalf_tilde[0]);
    T.insert(1, 3) = -0.25 * dt * dt / dx * hHalf_tilde[1];

    // All rows in between
    for (int i = 2; i < N; i++)
    {
        T.insert(i, i)
            = (epsilon * dx / (Gamma - 1.0)
               + (0.25 * dt * dt / dx)
                     * (hHalf_tilde[i] + hHalf_tilde[i - 1])); // diagonal component

        T.insert(i, i + 1)
            = -0.25 * dt * dt / dx * (hHalf_tilde[i] - hHalf_tilde[i - 1]); // upper 1
        T.insert(i, i + 2) = -0.25 * dt * dt / dx * hHalf_tilde[i];         // upper 2

        T.insert(i, i - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[i] - hHalf_tilde[i - 1]); // lower 1
        T.insert(i, i - 2) = -0.25 * dt * dt / dx * hHalf_tilde[i - 1];    // lower 2
    }

    // Row N-1
    T.insert(N - 1, N) = -0.25 * dt * dt / dx * hHalf_tilde[N - 1]
                         - 0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1]);
    T.insert(N - 1, N - 1)
        = epsilon * dx / (Gamma - 1.0)
          + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
    T.insert(N - 1, N - 2)
        = 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
    T.insert(N - 1, N - 3) = -0.25 * dt * dt / dx * hHalf_tilde[N - 2];

    // Row N
    T.insert(N, N) = -0.25 * dt * dt / dx * hHalf_tilde[N]
                     - 0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1])
                     + epsilon * dx / (Gamma - 1.0)
                     + 0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1]);
    T.insert(N, N - 1) = 0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1]);
    T.insert(N, N - 2) = -0.25 * dt * dt / dx * hHalf_tilde[N - 1];
}

void implicitUpdate(stateVec &F_imp, stateVec &consVec, stateVec &explicitConsVec,
                    double dt)
{
    int Niter = 2;
    double tol = 1e-6;

    std::vector<double> enthalpies(N);
    std::vector<double> hHalf_tilde(N);
    std::vector<double> initKineticEnergy(N);
    std::vector<double> initPressure(N);
    std::vector<double> initMomentum(N);

    initialisePicardVars(enthalpies, hHalf_tilde, initKineticEnergy, initPressure,
                         initMomentum, consVec, explicitConsVec);

    PicardSolver                     

    Eigen::SparseMatrix<double> T(N, N);
    setTMatrix(T, hHalf_tilde, dt);

    Eigen::VectorXd b(N);
    setbVector(b, consVec, explicitConsVec, hHalf_tilde, initKineticEnergy, dt);

    Eigen::GMRES<Eigen::SparseMatrix<double>> solver;
    solver.compute(T);
    solver.setMaxIterations(Niter);
    //solver.setTolerance(tol);

    Eigen::VectorXd pressure_sol = solver.solve(b);

    std::cout << "GMRES solved for pressure in " << Niter << " iterations" << std::endl;
}

void IMEXupdate(stateVec &consVec, stateVec &F_exp, stateVec &F_imp)
{

    double dt;
    compute_dt(consVec, dt);
    explicitUpdate(F_exp, consVec);
    stateVec explicitConsVec;
    implicitUpdate(F_imp, consVec, explicitConsVec, dt);

    for (int i = 1; i <= N; i++)
    {
        explicitConsVec[i][0] = F_exp[i][0];

        // consVec[i][1]
        //     = F_exp[i][1] + dt / (2 * epsilon * dx) * (p_iminus1 -
        //     p_iplus1);
    }
}
