#ifndef IMEX_H
#define IMEX_H

#include "parms.h"

#include "BoundaryConditions.h"
#include "EoS.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

typedef std::array<double, nVars> state;
typedef std::vector<state>        stateVec;

enum matrixBC
{
    transmissive,
    Dirichlet
};

void compute_dt(stateVec &consVec, double &dt, int iter);
void explicitFluxes(stateVec &consVec, stateVec &expFluxes);
void implicitFLuxes(stateVec &consVec, stateVec &impFluxes);
void explicitUpdate(stateVec &consVec, stateVec &expFluxes, double& dt);
void computeEnthalpyOnFace(std::vector<double> &hHalf_tilde,
                           std::vector<double> &enthalpies, stateVec &explicitConsVec);
void initialisePicardVars(std::vector<double> &enthalpies,
                          std::vector<double> &hHalf_tilde,
                          std::vector<double> &initKineticEnergy,
                          std::vector<double> &initPressure,
                          std::vector<double> &initMomentum, stateVec &consVec,
                          stateVec &explicitConsVec);
void setbVector(Eigen::VectorXd &b, stateVec &consVec, stateVec &explicitConsVec,
                std::vector<double> &hHalf_tilde, std::vector<double> &initKineticEnergy,
                double dt, matrixBC MatrixBC);
void setTMatrix(Eigen::SparseMatrix<double> &T, std::vector<double> &hHalf_tilde,
                double dt, matrixBC MatrixBC);
void computeImplicitPressure(std::vector<double> &implicitPressure, stateVec &consVec,
                             stateVec &explicitConsVec, double dt);
void computeImplicitEnthalpy(std::vector<double> &implicitEnthalpy,
                             std::vector<double>  implicitPressure);
void IMEXupdate(stateVec &consVec, double &dt);

#endif