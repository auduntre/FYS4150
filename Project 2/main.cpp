#include <iostream>

#include "jacobi.h"


int main(int argc, char **argv)
{
    arma::arma_rng::set_seed_random();

    uint n = 200;

    arma::mat A(n, n);
    arma::mat Ev(n, n);
    
    A.randu();
    A = A.t() * A;
    Ev.eye();

    arma::vec jacobival = jacobi(A, &Ev, 1.0E-12);

    return 0;
}