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
    arma::vec3 vel = {0, -2.0526e-6, 0};
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

    int matsize = (int) (0.05 * nTimesteps);
    arma::mat posMat(3, matsize, arma::fill::zeros); 
    arma::vec posNormMat(matsize, arma::fill::zeros);

    posMat.col(0) = mercury.position - sun.position;
    posNormMat(0) = arma::norm(posMat.col(0), 2);
    integrator.initialStep(sol);

    for (int i = 1; i <= nTimesteps - matsize; i++) {
        integrator.integrateOneStep(sol);
    }

    for (int i = 1; i < matsize; i++) {
        posMat.col(i) = mercury.position - sun.position;
        posNormMat(i) = arma::norm(posMat.col(i), 2);

        integrator.integrateOneStep(sol);
    }

    int nPeri = (int) (endTime / 0.24);
    int kPeri = 0;

    arma::mat periMat(3, nPeri, arma::fill::zeros);

    // Finding perihelion
    for (int i = 1; i < posNormMat.n_elem - 1; i++) {
        if (posNormMat(i) < posNormMat(i+1) && 
            posNormMat(i) < posNormMat(i-1)) 
        {
            if (kPeri < nPeri) {
                periMat.col(kPeri) = posMat.col(i);
                kPeri++;
            }
        }
    }

    periMat.save("../positions/mercury.txt", arma::raw_ascii);

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
