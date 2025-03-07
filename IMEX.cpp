#include "IMEX.h"
#include "BoundaryConditions.h"
#include "EoS.h"
#include "parms.h"
#include "unsupported/Eigen/IterativeSolvers"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale>
#include <set>
#include <vector>

extern double x0;
extern double x1;
extern double epsilon;
extern double Gamma;
extern double CFL;
extern int    N;
double        dx;

void compute_dt(stateVec &consVec, double &dt, int iter)
{

    double wave_speed;
    dx = (x1 - x0) / N;

    if (iter <= 3)
    {

        dx = (x1 - x0) / N;

        for (int i = 0; i < N; i++)
        {
            double rho      = consVec[i][0];
            double velocity = consVec[i][1] / rho;

            double max_speed
                = 2.0 * abs(velocity + sqrt(Gamma * getPressure(consVec[i]) / rho));

            if (max_speed > wave_speed)
            {
                wave_speed = max_speed;
            }
        }
    }

    else
    {
        for (int i = 0; i < N; i++)
        {
            double rho      = consVec[i][0];
            double velocity = consVec[i][1] / rho;

            double max_speed = 2.0 * abs(velocity);

            if (max_speed > wave_speed)
            {
                wave_speed = max_speed;
            }
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
        double v   = mom / rho;

        double k = getKineticEnergy(cons); // get kinetic energy

        expFlux[0] = mom;
        expFlux[1] = rho * v * v;
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

void explicitUpdate(stateVec &F_exp, stateVec &consVec, double &dt)
{

    stateVec expFluxes(N + nGhost);
    stateVec fPlusHalf(N + nGhost);
    stateVec fMinusHalf(N + nGhost);
    dx = (x1 - x0) / N;

    double aSpeed, aSpeedPlus, aSpeedMinus;

    explicitFluxes(consVec, expFluxes);

    for (int i = nGhost; i <= N; i++)
    {

        // aSpeed = std::max(std::fabs(consVec[i][1] / consVec[i][0]),
        //                   std::fabs(consVec[i + 1][1] / consVec[i + 1][0]));

        for (int var = 0; var < nVars; var++)
        {
            aSpeedPlus = std::max(fabs(consVec[i][1] / consVec[i][0]),
                                  fabs(consVec[i + 1][1] / consVec[i + 1][0]));
            aSpeedPlus = std::max(fabs(consVec[i][1] / consVec[i][0]),
                                  fabs(consVec[i - 1][1] / consVec[i - 1][0]));

            fPlusHalf[i][var]
                = 0.5 * (expFluxes[i + 1][var] + expFluxes[i][var])
                  - 0.5 * fabs(aSpeedPlus) * (consVec[i + 1][var] - consVec[i][var]);
            fMinusHalf[i][var]
                = 0.5 * (expFluxes[i][var] + expFluxes[i - 1][var])
                  - 0.5 * fabs(aSpeedMinus) * (consVec[i][var] - consVec[i - 1][var]);
        }
    }

    // Final explicit update
    for (int i = 2; i <= N; i++)
    {

        for (int var = 0; var < nVars; var++)
        {
            F_exp[i][var]
                = consVec[i][var] - dt / dx * (fPlusHalf[i][var] - fMinusHalf[i][var]);
        }
    }
}

void computeEnthalpyOnFace(std::vector<double> &hHalf_tilde,
                           std::vector<double> &enthalpies, stateVec &explicitConsVec)
{

    for (int i = 0; i < hHalf_tilde.size(); i++)
    {
        if (fabs(explicitConsVec[i][1]) <= 1e-10 && fabs(explicitConsVec[i + 1][1]) <= 1e-10)
        {
            hHalf_tilde[i] = fabs(0.5 * (enthalpies[i] + enthalpies[i + 1]));
        }

        else
        {

            hHalf_tilde[i] = fabs(0.5
                                  * ((enthalpies[i] * explicitConsVec[i][1])
                                     + (enthalpies[i + 1] * explicitConsVec[i + 1][1]))
                                  / (explicitConsVec[i][1] + explicitConsVec[i + 1][1]));
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

    computeEnthalpyOnFace(hHalf_tilde, enthalpies, explicitConsVec);

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
}

void setbVector(Eigen::VectorXd &b, stateVec &consVec, stateVec &explicitConsVec,
                std::vector<double> &hHalf_tilde, std::vector<double> &initKineticEnergy,
                double dt, matrixBC MatrixBC)
{

    // Elements outside the matrix bounds are considered to be the same as the elements at
    // the boundary Thus all coefficients are added together to the first element at the
    // boundary
    dx = (x1 - x0) / N;
    if (MatrixBC == transmissive)
    {
        b(0) = explicitConsVec[0][2] - initKineticEnergy[0]
               - 0.5 * dt / dx
                     * (hHalf_tilde[0] * explicitConsVec[1][1]
                        - hHalf_tilde[0] * explicitConsVec[0][1]);

        for (int i = 1; i < N - 1; i++)
        {
            b(i)
                = explicitConsVec[i][2] - initKineticEnergy[i]
                  - 0.5 * dt / dx
                        * (hHalf_tilde[i] * explicitConsVec[i + 1][1]
                           + (hHalf_tilde[i] - hHalf_tilde[i - 1]) * explicitConsVec[i][1]
                           - hHalf_tilde[i - 1] * explicitConsVec[i - 1][1]);
        }

        b(N - 1) = explicitConsVec[N - 1][2] - initKineticEnergy[N - 1]
                   - 0.5 * dt / dx
                         * (hHalf_tilde[N - 1] * explicitConsVec[N - 1][1]
                            - hHalf_tilde[N - 1 - 1] * explicitConsVec[N - 1 - 1][1]);

        b *= epsilon;
        b *= dx;
    }

    // The values outside the matrix bounds are set to a constant making the coefficient
    // constatns that can be move to the RHS of the equation - b
    else if (MatrixBC == Dirichlet)
    {
        double p_boundaryL = getPressure(consVec[0]);
        double p_boundaryR = getPressure(consVec[N - 1]);

        b(0) = explicitConsVec[0 + nGhost][2] - initKineticEnergy[0 + nGhost]
               - (0.5 * dt / dx
                  * (hHalf_tilde[0 + nGhost] * explicitConsVec[1 + nGhost][1]
                     - hHalf_tilde[0 + nGhost] * explicitConsVec[0 + nGhost][1]))
               - 0.25 * dt * dt / dx
                     * (hHalf_tilde[0 + nGhost] - hHalf_tilde[-1 + nGhost]) * p_boundaryL
               + 0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost] * p_boundaryL;

        b(1) = explicitConsVec[1 + nGhost][2] - initKineticEnergy[0 + nGhost]
               - (0.5 * dt / dx
                  * (hHalf_tilde[1 + nGhost] * explicitConsVec[1 + nGhost][1]
                     + explicitConsVec[1 + nGhost][1]
                           * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost])
                     - hHalf_tilde[0 + nGhost] * explicitConsVec[0 + nGhost][1]))
               + 0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost] * p_boundaryL;

        b(N - 1) = explicitConsVec[N - 1][2] - initKineticEnergy[N - 1]
                   - (0.5 * dt / dx
                      * (hHalf_tilde[N - 1] * explicitConsVec[N - 1][1]
                         + explicitConsVec[N - 1][1]
                               * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2])
                         - hHalf_tilde[N - 2] * explicitConsVec[N - 2][1]))
                   + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2])
                         * p_boundaryR
                   + 0.25 * dt * dt / dx * hHalf_tilde[N - 1]
                         * p_boundaryR; // out of bound so N = N-1

        b(N - 2) = explicitConsVec[N - 2][2] - initKineticEnergy[N - 2]
                   - (0.5 * dt / dx
                      * (hHalf_tilde[N - 2] * explicitConsVec[N - 2][1]
                         + (hHalf_tilde[N - 2] - hHalf_tilde[N - 3])
                               * explicitConsVec[N - 2][1]
                         - hHalf_tilde[N - 3] * explicitConsVec[N - 3][1]))
                   + 0.25 * dt * dt / dx * hHalf_tilde[N - 2] * p_boundaryR;

        b *= epsilon;
        b *= dx;
    }

    else
    {
        std::cout << "Invalid matrix boundary condition" << std::endl;
    }
    int j;
    for (int i = 2; i < N - 2; i++)
    {
        j    = i + nGhost;
        b(i) = explicitConsVec[j][2] - initKineticEnergy[j]
               - (0.5 * dt / dx
                  * (hHalf_tilde[j] * explicitConsVec[j + 1][1]
                     + (hHalf_tilde[j] - hHalf_tilde[j - 1]) * explicitConsVec[j][1]
                     - hHalf_tilde[j - 1] * explicitConsVec[j - 1][1]));
    }
}

void setTMatrix(Eigen::SparseMatrix<double> &T, std::vector<double> &hHalf_tilde,
                double dt, matrixBC MatrixBC)
{

    dx = (x1 - x0) / N;
    // setting up matrix T using eigen
    T.reserve(Eigen::VectorXi::Constant(N, 5));

    // Manually account for boundary conditions for fist and last two rows
    // Transmissive case

    if (MatrixBC == transmissive)
    {
        // Row 1
        T.insert(0, 0)
            = -0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost]
              + 0.25 * dt * dt / dx * (hHalf_tilde[0 + nGhost] - hHalf_tilde[0 + nGhost])
              + epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[0 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(0, 1)
            = -0.25 * dt * dt / dx * (hHalf_tilde[0 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(0, 2) = -0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost];

        // Row 2
        T.insert(1, 0)
            = -0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost]
              - 0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 1)
            = epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 2)
            = -0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 3) = -0.25 * dt * dt / dx * hHalf_tilde[1 + nGhost];

        // Row N-1
        T.insert(N - 1 - 1, N - 1)
            = -0.25 * dt * dt / dx * hHalf_tilde[N - 1]
              - 0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1]);
        T.insert(N - 1 - 1, N - 1 - 1)
            = epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
        T.insert(N - 1 - 1, N - 2 - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
        T.insert(N - 1 - 1, N - 3 - 1) = -0.25 * dt * dt / dx * hHalf_tilde[N - 2];

        // Row N
        T.insert(N - 1, N - 1)
            = -0.25 * dt * dt / dx * hHalf_tilde[N - 1]
              - 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 1 - 1])
              + epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 1 - 1]);
        T.insert(N - 1, N - 1 - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 1 - 1]);
        T.insert(N - 1, N - 2 - 1) = -0.25 * dt * dt / dx * hHalf_tilde[N - 1 - 1];
    }

    else if (MatrixBC == Dirichlet)
    {
        // Row 1
        T.insert(0, 0)
            = epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[0 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(0, 1)
            = -0.25 * dt * dt / dx * (hHalf_tilde[0 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(0, 2) = -0.25 * dt * dt / dx * hHalf_tilde[0 + nGhost];

        // Row 2
        T.insert(1, 0)
            = -0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 1)
            = epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 2)
            = -0.25 * dt * dt / dx * (hHalf_tilde[1 + nGhost] - hHalf_tilde[0 + nGhost]);
        T.insert(1, 3) = -0.25 * dt * dt / dx * hHalf_tilde[1 + nGhost];

        // Row N-1
        T.insert(N - 1 - 1, N - 1)
            = -0.25 * dt * dt / dx * (hHalf_tilde[N] - hHalf_tilde[N - 1]);
        T.insert(N - 1 - 1, N - 1 - 1)
            = epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
        T.insert(N - 1 - 1, N - 2 - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 2]);
        T.insert(N - 1 - 1, N - 3 - 1) = -0.25 * dt * dt / dx * hHalf_tilde[N - 2];

        // Row N
        T.insert(N - 1, N - 1)
            = +epsilon * dx / (Gamma - 1.0)
              + 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 1 - 1]);
        T.insert(N - 1, N - 1 - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[N - 1] - hHalf_tilde[N - 1 - 1]);
        T.insert(N - 1, N - 2 - 1) = -0.25 * dt * dt / dx * hHalf_tilde[N - 1 - 1];
    }

    else
    {
        std::cout << "Invalid matrix boundary condition" << std::endl;
    }

    // All rows in between
    int j;
    for (int i = 2; i < N - 2; i++)
    {
        j = i + nGhost;

        T.insert(i, i)
            = (epsilon * dx / (Gamma - 1.0)
               + (0.25 * dt * dt / dx)
                     * (hHalf_tilde[j] + hHalf_tilde[j - 1])); // diagonal component

        T.insert(i, i + 1)
            = -0.25 * dt * dt / dx * (hHalf_tilde[j] - hHalf_tilde[j - 1]); // upper 1
        T.insert(i, i + 2) = -0.25 * dt * dt / dx * hHalf_tilde[j];         // upper 2

        T.insert(i, i - 1)
            = 0.25 * dt * dt / dx * (hHalf_tilde[j] - hHalf_tilde[j - 1]); // lower 1
        T.insert(i, i - 2) = -0.25 * dt * dt / dx * hHalf_tilde[j - 1];    // lower 2
    }
}

void computeImplicitPressure(std::vector<double> &implicitPressure, stateVec &consVec,
                             stateVec &explicitConsVec, double dt)
{
    int      Niter    = 2;
    double   tol      = 1e-6;
    matrixBC MatrixBC = Dirichlet;

    std::vector<double> enthalpies(N + nGhost);
    std::vector<double> hHalf_tilde(N + nGhost);
    std::vector<double> initKineticEnergy(N + nGhost);
    std::vector<double> initPressure(N + nGhost);
    std::vector<double> initMomentum(N + nGhost);

    initialisePicardVars(enthalpies, hHalf_tilde, initKineticEnergy, initPressure,
                         initMomentum, consVec, explicitConsVec);

    Eigen::SparseMatrix<double> T(N, N);
    setTMatrix(T, hHalf_tilde, dt, MatrixBC);

    // std::cout << Eigen::MatrixXd(T) << std::endl;

    Eigen::VectorXd b(N);
    setbVector(b, consVec, explicitConsVec, hHalf_tilde, initKineticEnergy, dt, MatrixBC);

    // std::cout << Eigen::VectorXd(b) << std::endl;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    // Eigen::GMRES<Eigen::SparseMatrix<double> > solver;
    solver.compute(T);
    solver.setMaxIterations(Niter);
    // solver.setTolerance(tol);

    Eigen::VectorXd pressure_sol = solver.solve(b);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "GMRES failed to converge" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = nGhost; i < N; i++)
    {
        implicitPressure[i] = pressure_sol(i);
    }

    transmissiveBC(implicitPressure);
}

void computeImplicitEnthalpy(std::vector<double> &implicitEnthalpy,
                             std::vector<double>  implicitPressure)
{

    for (int i = 0; i < N; i++)
    {
        implicitEnthalpy[i] = implicitPressure[i] * Gamma / (Gamma - 1.0);
    }
}

void IMEXupdate(stateVec &consVec, double &dt)
{
    std::vector<double> implicitPressure(N + nGhost);
    std::vector<double> implicitEnthalpy(N + nGhost);
    std::vector<double> implicit_hTilde(N + nGhost);
    stateVec            F_exp(N + nGhost), F_imp(N + nGhost);
    dx = (x1 - x0) / N;

    // Do explicit update fist
    explicitUpdate(F_exp, consVec, dt);
    transmissiveBC(F_exp);

    // get implicit variables
    computeImplicitPressure(implicitPressure, consVec, F_exp, dt);
    computeImplicitEnthalpy(implicitEnthalpy, implicitPressure);

    for (int i = nGhost; i <= N; i++)
    {
        consVec[i][0] = F_exp[i][0];

        consVec[i][1] = F_exp[i][1]
                        - 0.5 * dt / (epsilon * dx)
                              * (implicitPressure[i + 1] - implicitPressure[i - 1]);

        computeEnthalpyOnFace(implicit_hTilde, implicitEnthalpy, F_exp); 
    }

    transmissiveBC(consVec);

    for (int i = nGhost; i <= N; i++)
    {
        consVec[i][2]
            = F_exp[i][2]
              - 0.5 * dt / dx
                    * (implicit_hTilde[i] * (consVec[i + 1][1] + consVec[i][1])
                       - implicit_hTilde[i - 1] * (consVec[i][1] + consVec[i - 1][1]));
    }
}
