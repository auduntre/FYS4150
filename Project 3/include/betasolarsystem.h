#ifndef BETASOLARSYSTEM_H
#define BETASOLARSYSTEM_H

#include "solarsystem.h"


class BetaSolarSystem : public SolarSystem
{
    private:
        double beta;

    public:
        BetaSolarSystem (double b);

        double forcePower (double r, CelestialBody body);
}; 

#endif
