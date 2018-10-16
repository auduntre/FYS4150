#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>

#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"

int main (int argc, char **argv)
{
    void bodyPrint (std::vector<CelestialBody> &bodies);

    double endTime = 1.00;
    int nTimesteps = 1000;

    if (argc >= 2) {
        nTimesteps = std::atoi(argv[1]);

        if (argc >= 3) {
            endTime = std::atof(argv[2]);
        }
    }

    SolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, 0, 0};
    double mass = 1.0;

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {0.3075, 0, 0};
    vel = {0, 12.44, 0};
    mass  = 1.65e-7;

    CelestialBody &mercury = sol.createCelestialBody(pos, vel, mass);

    std::vector<CelestialBody> &bodies = sol.bodies();
    bodyPrint(bodies);


    double dt = endTime / nTimesteps;

    std::cout << "################# CALCULATE #################" << std::endl;
    Verlet integrator(dt);


    arma::mat posMat(3, nTimesteps, arma::fill::zeros); 
    posMat.col(0) = mercury.position;
    integrator.initialStep(sol);

    for (int i = 1; i < nTimesteps; i++) {
        posMat.col(i) = mercury.position;
        integrator.integrateOneStep(sol);
    }

    posMat.save("../positions/mercury.txt", arma::raw_ascii);

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
