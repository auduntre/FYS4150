#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>

#include "jacobi.h"
#include "toeplitz.h"


int main(int argc, char **argv)
{
    arma::vec analy_eigval(uint N);
    
    uint N = 100;
    
    if (argc > 1) {
        N = std::atoi(argv[1]);
    }

    // Value thats spits out good approx to analytical eigenvalues
    double rhoN = 5.00;

    if (argc > 2) {
        rhoN = std::atof(argv[2]);
    }

    double h = rhoN  / N;

    double dconst = 2.0 / (h * h);
    double econst = -1.0 / (h * h);

    double eps = 1E-8;

    arma::vec rho = arma::linspace<arma::vec>(h, rhoN - h, N-1);
    arma::vec Vi = (rho % rho);
    
    arma::vec eigval = analy_eigval(N);
    arma::vec jeigval(N-1);
    arma::vec tmp(4);

    arma::mat A(N-1, N-1);
    arma::mat Ev = arma::eye<arma::mat>(N-1, N-1);

    std::string tmpname= "results/approxN" + std::to_string(N) + \
                         "rhoN" + std::to_string((uint) rhoN) + ".txt";

    Tridiag td(econst, dconst, N-1);
    td.set_variable_diags(Vi);
    A = td.get_matrix();

    jeigval = arma::sort(jacobi(A, &Ev, eps));

    for (uint i = 0; i < tmp.n_elem; i++) {
        tmp(i) = std::fabs(jeigval(i) - eigval(i));
    }

    tmp.save(tmpname, arma::raw_ascii);

    return 0;
}


arma::vec analy_eigval(uint N)
{
    return arma::linspace<arma::vec>(3.0, 3.0 + 4.0*(N-2), N-1);
}