#include <iostream>
#include <cstdlib>
#include <iomanip> 

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
    arma::vec3 vel = {0, 0, 0};
    double mass = 1.0;

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {1, 0, 0};
    vel = {0, 2*M_PI, 0};
    mass  = 1e-6;

    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    std::vector<CelestialBody> &bodies = sol.bodies();
    bodyPrint(bodies);


    double endTime = 1.0;
    double dt = endTime / nTimesteps;

    std::cout << "################# CALCULATE #################" << std::endl;
    Verlet integrator(dt);
    integrator.integrateNtimes(sol, nTimesteps);

    bodyPrint(bodies);
    
    sol.writeToFile("positions.xyz");

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