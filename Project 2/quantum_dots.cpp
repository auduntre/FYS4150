#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <string>

#include "jacobi.h"
#include "toeplitz.h"
#include "haroscill.h"


int main (int argc, char **argv)
{
    arma::vec one_electron (uint N, double rhoN);
    void two_electrons (uint N, double rhoN);

    double omegas[4] = {0.01, 0.50, 1.00, 5.00};
    uint ns[7] = {10, 25, 50, 75, 100, 150, 200};

    uint N = 100;
    uint jmin;

    // Value thats spits out good approx to analytical eigenvalues
    double rhoN = 5.00;

    if (argc > 1) {
        N = std::atoi(argv[1]);
    }

    if (argc > 2) {
        rhoN = std::atof(argv[2]);
    }
    
    arma::mat Ev(N-1, N-1, arma::fill::eye);
    arma::vec jeigval(N-1);
    arma::vec tmp;
    arma::vec evmin;

    std::string tmpname;
    std::string eiglname; // eigenvalue
    std::string eigcname; // eigenvector
    std::ofstream eigvaluefile;

    // For one electron how close is the eigenvalues for diffirent Ns
    for (uint n: ns) {
        Harmonic_oscillator hon = Harmonic_oscillator(rhoN, n);

        tmpname = "results/approxN" + std::to_string(n) + \
                  "rhoN" + std::to_string((uint) rhoN) + ".txt";

        tmp = hon.one_electron();
        tmp.save(tmpname, arma::raw_ascii);
    }

    Harmonic_oscillator ho = Harmonic_oscillator(rhoN, N);

    // Run two electron case for diffirent omega values
    for (double omega: omegas) {
        Ev.eye();

        eigcname = "results/eigvector" + std::to_string(N) + \
                   "rhoN" + std::to_string((uint) rhoN) + \
                   "omega" + std::to_string(omega) + ".txt";

        eiglname = "results/eigvalue" + std::to_string(N) + \
                   "rhoN" + std::to_string((uint) rhoN) + \
                   "omega" + std::to_string(omega) + ".txt";

        jeigval = ho.two_electrons(&Ev, omega, true);
        jmin = jeigval.index_min();

        eigvaluefile.open(eiglname);
        eigvaluefile << jeigval(jmin) << std::endl;
        eigvaluefile.close();

        evmin = Ev.col(jmin);
        evmin.save(eigcname, arma::raw_ascii);
    }

    Ev.eye();

    // Run case with non-interacting system 
    jeigval = ho.two_electrons(&Ev, omegas[0], false);
    jmin = jeigval.index_min();

    eigvaluefile.open("results/noninteracting.txt");
    eigvaluefile << jeigval(jmin) << std::endl;
    eigvaluefile.close();

    evmin = Ev.col(jmin);
    evmin.save("results/noninteractionvector.txt", arma::raw_ascii);

    return 0;
}
