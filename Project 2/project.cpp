#include <armadillo>
#include <cmath>

#include "project.h"

arma::uword maxoffdiag (arma::mat X)
{
    // Set diag to zero and find index of extreme value in matrix
    return arma::abs(arma::trimatu(X)).index_max();
}
