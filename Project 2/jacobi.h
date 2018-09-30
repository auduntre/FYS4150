#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>

typedef unsigned uint;

void maxoffdiag (arma::mat X);

arma::mat jacobi (arma::mat A, double eps, uint maxiter);

#endif
