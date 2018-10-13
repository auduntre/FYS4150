#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.h"

#include <fstream>
#include <vector>
#include <string>
#include <cmath>

//static const double GravConst = 6.67408e-11;
static const double GravConst = 4.0 * M_PI * M_PI; // Scaled


class SolarSystem
{
    private:
        std::vector<CelestialBody> bods;
        std::ofstream ofile;
        
        arma::vec3 angMom;
        double ke;  // Kinetic Energy
        double pe;  // Potential Energy
        double beta; // force-beta

    public:
        SolarSystem();

        CelestialBody &createCelestialBody (arma::vec3 position, arma::vec3 velocity,
                                            double mass);
        
        void calculateForcesAndEnergy ();
        int numberOfBodies () const;
        
        double totalEnergy () const;
        double potentialEnergy () const;
        double kineticEnergy () const;
        
        void writeToFile (std::string filename);
        
        arma::vec3 angularMomentum () const;
        std::vector<CelestialBody> &bodies ();
};

#endif
