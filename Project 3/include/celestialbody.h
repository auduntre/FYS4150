#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include <armadillo>


class CelestialBody
{
    public:
        arma::vec3 position;
        arma::vec3 velocity;
        arma::vec3 force;
        arma::vec3 orbAngMom;
        
        double mass;
        
        arma::vec3 past_force;

        CelestialBody (arma::vec3 pos, arma::vec3 vel, double m);
        CelestialBody (double x, double y, double z, double vx, double vy,
                       double vz, double mass);

        void computeOrbitalAngMom (); 
        void resetForce ();
};

#endif
