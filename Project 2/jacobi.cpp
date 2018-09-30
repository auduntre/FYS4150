#include <iostream>
#include <cmath>

#include "jacobi.h"


arma::mat jacobi (arma::mat A, double eps, uint maxiter) 
{
    double maxoffdiag (arma::mat X, uint *i, uint *j);
    arma::mat rotate (arma::mat B, uint k, uint l);

    arma::mat B = A;

    double maxval;

    uint maxi;
    uint maxj;
    uint iter = 0;

    maxval = maxoffdiag(B, &maxi, &maxj);

    while (maxval > eps && maxiter > iter) {
        B = rotate(B, maxi, maxj);
        maxval = maxoffdiag(B, &maxi, &maxj);
        iter++;
    }

    return B;
}


arma::mat rotate (arma::mat B, uint k, uint l)
{
    arma::vec aik = B.col(k);
    arma::vec ail = B.col(l);
    arma::uvec indices;

    double akl = B(k, l);
    double akk = B(k, k);
    double all = B(l, l); 

    double sign;
    double tau;
    double tanw;
    double sinw;
    double cosw;

    indices = arma::linspace<arma::uvec>(0, B.n_rows - 1, B.n_rows);
    
    //Not iterate thru these indices
    indices.shed_row(k);
    indices.shed_row(l - 1);

    if (akl != 0.0 ) {   
        tau = (all - akk) / (2 * akl);
        sign = tau / std::fabs(tau);

        tanw = sign / (std::fabs(tau) + std::sqrt(1.0 + tau*tau));

        cosw = 1 / std::sqrt(1 + tanw * tanw);
        sinw = cosw*tanw;
    } else {
        cosw = 1.0;
        sinw = 0.0;
    }

    B(k,k) = cosw * cosw * akk - 2.0 * cosw * sinw * akl + sinw * sinw * all;
    B(l,l) = sinw * sinw * akk + 2.0 * cosw * sinw * akl + cosw * cosw * all;
    
    // hard-coding non-diagonal elements
    B(k,l) = 0.0;
    B(l,k) = 0.0;
    
    for (arma::uword i: indices) {
            B(i,k) = cosw * aik(i) - sinw * ail(i);
            B(k,i) = B(i,k);
            B(i,l) = cosw * ail(i) + sinw * aik(i);
            B(l,i) = B(i,l);
    }

    return B;
}


double maxoffdiag (arma::mat X, uint *i, uint *j)
{
    // Set diag to zero and find index of extreme value in matrix
    uint ij = arma::abs(arma::trimatu(X) - arma::diagmat(X)).index_max();

    *j = ij / X.n_rows;
    *i = ij - (*j) * X.n_rows;

    return std::fabs(X(ij));
}