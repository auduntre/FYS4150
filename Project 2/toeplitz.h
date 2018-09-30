#ifndef TOEPLITZ_H
#define TOEPLITZ_H

#include <armadillo>

typedef unsigned uint;

arma::mat sym_tridiag (double a, double d, uint N);

#endif