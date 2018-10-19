#include "gensolarsystem.h"

GenSolarSystem::GenSolarSystem ()
{
}


double GenSolarSystem::forcePower (double r, CelestialBody body)
{
    double l = arma::norm(body.orbAngMom, 2);
    return r / (1.0 + 3.0 * l * l / (r * r * c * c)); 
}
