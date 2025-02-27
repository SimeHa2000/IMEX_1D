#ifndef PARMS_H
#define PARMS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

const int nVars = 3; // Define nVars here

typedef std::array<double, nVars> state;
typedef std::vector<state> stateVec;

extern double t_stop;
extern int N;
extern double Gamma;
extern double CFL;
extern double epsilon;
extern double x0;
extern double x1;
extern int nGhost;

std::map<std::string, std::string> readSettingsFile(std::string filename);
int readParams(int argc,char** argv);
void readStates(state &stateL, state &stateR, int argc,char** argv);

#endif