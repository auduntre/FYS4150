#include <ctime>
#include <iostream>

#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"


int main(int argc, char **argv)
{
    SolarSystem createSystem ();

    int nTimeSteps = 100000000;
    double endTime = 10;

    double dt = endTime / nTimeSteps;

    clock_t start;
    clock_t end;
    
    Euler integratorE(dt);
    Verlet integratorV(dt);
    
    SolarSystem solV = createSystem();

    // MEASURE TIME VERLET
    start = std::clock();
    integratorV.integrateNtimes(solV, nTimeSteps);
    end = std::clock();

    std::cout << "Time for Verlet: " << (double) (end - start) / CLOCKS_PER_SEC
              << std::endl;

    SolarSystem solE = createSystem();
    
    // MEASURE TIME EULER
    start = std::clock();
    integratorE.integrateNtimes(solE, nTimeSteps);
    end = std::clock();

    std::cout << "Time for Euler: " << (double) (end - start) / CLOCKS_PER_SEC
              << std::endl;
    
    return 0;
}


SolarSystem createSystem ()
{
    SolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, -1.8849555921538758e-05, 0};
    double mass = 1.0;
    
    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);
    
    pos = {1, 0, 0};
    vel = {0, 8.9, 0};
    mass  = 3e-6;
    
    CelestialBody &earth = sol.createCelestialBody(pos, vel, mass);

    return sol;
}
