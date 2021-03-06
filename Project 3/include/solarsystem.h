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

    public:
        SolarSystem();

        CelestialBody &createCelestialBody (arma::vec3 position, arma::vec3 velocity,
                                            double mass);
        
        void calculateForcesAndEnergy (bool forceAndPE, bool angMomAndKE);
        virtual double forcePower (double r, CelestialBody body);

        // Getters
        int numberOfBodies () const;
        double totalEnergy () const;
        double potentialEnergy () const;
        double kineticEnergy () const;
        arma::vec3 angularMomentum () const;
        std::vector<CelestialBody> &bodies ();
        
        void writeToFile (std::string filename);
        
};

#endif
