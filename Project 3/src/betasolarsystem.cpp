#include <cmath>

#include "betasolarsystem.h"

BetaSolarSystem::BetaSolarSystem (double b) : beta(b - 1.0)
{
}


double BetaSolarSystem::forcePower (double r)
{
    return std::pow(r, this->beta);
}
