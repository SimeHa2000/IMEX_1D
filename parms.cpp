#include "parms.h"

double t_stop;
int N;
double Gamma;
double CFL;
double epsilon;
double x0;
double x1;  
int nGhost;

std::map<std::string, std::string> readSettingsFile(std::string filename){
    std::map<std::string, std::string> settings;
    std::ifstream settingsFile(filename);
    
    if (settingsFile.is_open()){
        std::string line;

        while (std::getline(settingsFile, line)){

            if(line.empty() || line[0] == '#'){
                continue;
            }

            std::istringstream iss(line);
            std::string key;

            if (std::getline(iss, key, '=')) {
                std::string value;
                if (std::getline(iss, value)) {

                    // Remove spaces around key and value
                    key.erase(0, key.find_first_not_of(" \t"));
                    key.erase(key.find_last_not_of(" \t") + 1);
                    value.erase(0, value.find_first_not_of(" \t"));
                    value.erase(value.find_last_not_of(" \t") + 1);
                    
                    settings[key] = value;
                }
            }
        }
    }

    else{
        std::cout << "Unable to open settings file" << std::endl;
    }

    return settings;
}

int readParams(int argc,char** argv){

    if(argc < 2){
        std::cout << "Please provide a settings file" << std::endl;
        return 1;
    }

    std::map<std::string, std::string> settings = readSettingsFile(argv[1]);
    
    t_stop = std::stod(settings["t_stop"]); 
    N = std::stoi(settings["N"]);
    Gamma = std::stod(settings["Gamma"]);
    CFL = std::stod(settings["CFL"]);
    epsilon = std::stod(settings["epsilon"]);
    x0 = std::stod(settings["x0"]);
    x1 = std::stod(settings["x1"]);
    nGhost = std::stoi(settings["nGhost"]);
    

    return 0;
}

void readStates(state &stateL, state &stateR, int argc,char** argv){

    std::map<std::string, std::string> settings = readSettingsFile(argv[1]);

    double rhoL = std::stod(settings["rhoL"]);
    double vL = std::stod(settings["vL"]);
    double pL = std::stod(settings["pL"]);

    double rhoR = std::stod(settings["rhoR"]);
    double vR = std::stod(settings["vR"]);
    double pR = std::stod(settings["pR"]);

    stateL[0] = rhoL;
    stateL[1] = vL;
    stateL[2] = pL;

    stateR[0] = rhoR;   
    stateR[1] = vR;
    stateR[2] = pR;
}