#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

#include "jacobi.h"
#include "toeplitz.h"
#include "haroscill.h"


int main (int argc, char **argv)
{
    arma::vec one_electron (uint N, double rhoN);
    void two_electrons (uint N, double rhoN);

    uint N = 100;
    double rhoN = 5.00;
    
    arma::vec tmp;
    std::string tmpname;

    if (argc > 1) {
        N = std::atoi(argv[1]);
    }

    // Value thats spits out good approx to analytical eigenvalues
    
    if (argc > 2) {
        rhoN = std::atof(argv[2]);
    }

    Harmonic_oscillator ho = Harmonic_oscillator(rhoN, N);

    tmpname = "results/approxN" + std::to_string(N) + \
              "rhoN" + std::to_string((uint) rhoN) + ".txt";


    // = ho.one_electron();
    tmp.save(tmpname, arma::raw_ascii);

    return 0;
}
