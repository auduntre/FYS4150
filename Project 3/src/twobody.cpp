#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <string>

#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"

int main (int argc, char **argv)
{
    void bodyPrint (std::vector<CelestialBody> &bodies);

    int nTimesteps = 1000;

    if (argc >= 2) {
        nTimesteps = std::atoi(argv[1]);
    }

    SolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, -1.8849555921538758e-05, 0};
    double mass = 1.0;

    if (argc >= 4) {
        vel(1) = - std::atof(argv[3]) * 3e-6;
    }
    
    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {1, 0, 0};
    vel = {0, 2 * M_PI, 0};
    mass  = 3e-6;

    if (argc >= 4) {
        vel(1) = std::atof(argv[3]);
    }
    
    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    std::vector<CelestialBody> &bodies = sol.bodies();
    bodyPrint(bodies);


    double endTime = 10.0;

    if (argc >= 3) {
        endTime = std::atof(argv[2]);
    }
    double dt = endTime / nTimesteps;

    std::cout << "################# CALCULATE #################" << std::endl;
    Verlet integrator(dt);
    integrator.initialStep(sol);
    
    arma::mat posMat(3, nTimesteps, arma::fill::zeros); 
    posMat.col(0) = earth.position;

    for (int i = 1; i < nTimesteps; i++) {
        integrator.integrateOneStep(sol);

        posMat.col(i) = earth.position;
    }

    if (argc >= 4) {
        std::string velstr(argv[3]);
        posMat.save("../positions/escape" + velstr + ".txt",
                     arma::raw_ascii);
    } else {
        posMat.save("../positions/positions" + std::to_string(nTimesteps) + ".txt",
                    arma::raw_ascii);
    }

    bodyPrint(bodies);
    
    return 0;
}


void bodyPrint (std::vector<CelestialBody> &bodies)
{
    for(int i = 0; i < bodies.size(); i++) {
        CelestialBody &body = bodies[i];

        std::cout << "--------------------------------" 
                  << "--------------------------------" << std::endl;

        std::cout << "The position of this object is: " << body.position.t()
                  << "The velocity of this object is: " << body.velocity.t()
                  << std::endl;
    }
}
