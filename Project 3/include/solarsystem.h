#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.h"

#include <fstream>
#include <vector>
#include <string>


static const double GravConst = 6.67408e-11;


class SolarSystem
{
    private:
        std::vector<CelestialBody> bods;
        std::ofstream ofile;
        
        vec3 angMom;
        double ke;  // Kinetic Energy
        double pe;  // Potential Energy

    public:
        SolarSystem();

        CelestialBody &createCelestialBody (vec3 position, vec3 velocity,
                                            double mass);
        
        void calculateForcesAndEnergy ();
        int numberOfBodies () const;
        
        double totalEnergy () const;
        double potentialEnergy () const;
        double kineticEnergy () const;
        
        void writeToFile (std::string filename);
        
        arma::vec angularMomentum () const;
        std::vector<CelestialBody> &bodies ();
};

#endif
