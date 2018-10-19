#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>

#include "solarsystem.h"
#include "gensolarsystem.h"
#include "euler.h"
#include "verlet.h"

int main (int argc, char **argv)
{
    void bodyPrint (SolarSystem &sol);

    double endTime = 1.00;
    int nTimesteps = 1000;

    if (argc >= 2) {
        nTimesteps = std::atoi(argv[1]);

        if (argc >= 3) {
            endTime = std::atof(argv[2]);
        }
    }

    GenSolarSystem sol;

    arma::vec3 pos = {0, 0, 0};
    arma::vec3 vel = {0, -2.0651820597435255e-06, 0};
    double mass = 1.0;

    CelestialBody &sun = sol.createCelestialBody(pos, vel, mass);

    pos = {0.3075, 0, 0};
    vel = {0, 12.44, 0};
    mass  = 1.66011419593531e-07;

    CelestialBody &mercury = sol.createCelestialBody(pos, vel, mass);

    double dt = endTime / nTimesteps;

    bodyPrint(sol);
    std::cout << "################# CALCULATE #################" << std::endl;

    Verlet integrator(dt);
    integrator.initialStep(sol);
    bodyPrint(sol);

    int nPeri = (int) (endTime / 0.24);
    int kPeri = 0;

    arma::mat periMat(3, nPeri, arma::fill::zeros);
    
    arma::vec3 futurePos = mercury.position - sun.position;
    arma::vec3 nowPos = futurePos  + 1;
    
    double futureNorm = arma::norm(futurePos, 2);
    double nowNorm = futureNorm + 1;
    double pastNorm = futureNorm + 1;

    for (int i = 1; i <= nTimesteps; i++) {
        integrator.integrateOneStep(sol);

        nowPos = futurePos;
        futurePos = mercury.position - sun.position;

        pastNorm = nowNorm;
        nowNorm = futureNorm;
        futureNorm = arma::norm(nowPos, 2);

        if (nowNorm < futureNorm && nowNorm < pastNorm) {
            periMat.col(kPeri) = nowPos;
            kPeri++;
        }
    }

    periMat.save("../positions/mercury.txt", arma::raw_ascii);
    bodyPrint(sol);

    return 0;
}


void bodyPrint (SolarSystem &sol)
{
    std::vector<CelestialBody> &bodies = sol.bodies();

    std::cout << "The Systems AngMom:" << sol.angularMomentum().t();
    std::cout << "The Systems KE:    " << sol.kineticEnergy() << std::endl;
    std::cout << "The Systems PE:    " << sol.potentialEnergy() << std::endl;
    std::cout << "The Systems TE:    " << sol.totalEnergy() << std::endl;

    for(int i = 0; i < bodies.size(); i++) {
        CelestialBody &body = bodies[i];

        std::cout << "--------------------------------" 
                  << "--------------------------------" << std::endl;

        std::cout << "The position of this object is: " << body.position.t()
                  << "The velocity of this object is: " << body.velocity.t()
                  << std::endl;
    }
}
