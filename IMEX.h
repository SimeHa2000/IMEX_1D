#ifndef IMEX_H
#define IMEX_H

#include "parms.h"

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Spectra/SymEigsSolver.h>


typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

void compute_dt(stateVec& consVec, double& dt);
void explicitFluxes(stateVec& consVec, stateVec& expFluxes);
void implicitFLuxes(stateVec& consVec, stateVec& impFluxes);
void explicitUpdate(stateVec& consVec, stateVec& expFluxes, double dt);
void initialisePicardVars(std::vector<double>& enthalpies, std::vector<double>& hHalf_tilde, std::vector<double>& initKineticEnergy, std::vector<double>& initPressure, std::vector<double>& initMomentum, stateVec& consVec, stateVec& explicitConsVec);
void setbVector(Eigen::VectorXd& b, stateVec& consVec, stateVec& explicitConsVec, std::vector<double>& hHalf_tilde, std::vector<double>& initKineticEnergy, double dt);
void setTMatrix(Eigen::SparseMatrix<double>& T, std::vector<double>& hHalf_tilde, double dt);
void computeImplicitPressure(std::vector<double>& implicitPressure, stateVec& consVec, stateVec& explicitConsVec, stateVec& F_exp, double dt);
void IMEXupdate(stateVec& consVec, stateVec& F_exp, stateVec& F_imp);

#endif