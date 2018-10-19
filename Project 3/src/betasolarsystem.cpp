#include <cmath>

#include "betasolarsystem.h"

BetaSolarSystem::BetaSolarSystem (double b) : beta(b - 1.0)
{
}


double BetaSolarSystem::forcePower (double r, CelestialBody body)
{
    return std::pow(r, this->beta);
}
