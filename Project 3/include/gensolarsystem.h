#ifndef GENSOLARSYSTEM_H
#define GENSOLARSTSTEM_H

#include "solarsystem.h"

// Speed of light in AU / years
static const double c = 63239.7263

class GenSolarSystem : public SolarSystem
{
    private:
        double l;

    public:
        GenSolarSystem ();

        double forcePower (double r);
}

#endif
