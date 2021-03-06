#include <iostream>

#include "jacobi.h"


arma::vec jacobi (arma::mat A, arma::mat *Ev, double eps) 
{
    double maxoffdiag (arma::mat X, uint *i, uint *j);
    arma::mat rotate (arma::mat A, arma::mat *Ev, uint k, uint l);

    double maxval;

    uint maxi;
    uint maxj;
    uint maxiter = A.n_rows * A.n_rows * A.n_rows;
    uint iter = 0;

    maxval = maxoffdiag(A, &maxi, &maxj);
    
    // The jacobi's rotation algorithm
    while (maxval > eps && maxiter > iter) {
        A = rotate(A, Ev, maxi, maxj);
        maxval = maxoffdiag(A, &maxi, &maxj);
        iter++;
    }

    // Give number of iterations and wheter sucsess
    if (iter >= maxiter) {
        std::cout << "Did not converge to max offdiagnonal element <= " << eps
                  << " for " << maxiter << " iterations." << std::endl;
    } else {
        std::cout << "Solution found in " << iter 
                  << " iterations." << std::endl;
    }
    
    return A.diag();
}


arma::mat rotate (arma::mat B, arma::mat *Ev, uint k, uint l)
{
    arma::mat R = *Ev;

    arma::vec aik = B.col(k);
    arma::vec ail = B.col(l);
    arma::vec rik = R.col(k);
    arma::vec ril = R.col(l);
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

    // Rotational elements
    if (akl != 0.0 ) {   
        tau = (all - akk) / (2 * akl);
        
        if (tau > 0.0) {
            tanw = 1.0 / (tau + std::sqrt(1.0 + tau*tau));
        } else {
            tanw = -1.0 / (-tau + std::sqrt(1.0 + tau*tau));
        }

        cosw = 1 / std::sqrt(1 + tanw * tanw);
        sinw = cosw*tanw;
    } else {
        cosw = 1.0;
        sinw = 0.0;
    }

    // Values on the diagonal
    B(k,k) = cosw * cosw * akk - 2.0 * cosw * sinw * akl + sinw * sinw * all;
    B(l,l) = sinw * sinw * akk + 2.0 * cosw * sinw * akl + cosw * cosw * all;
    
    // hard-coding non-diagonal elements
    B(k,l) = 0.0;
    B(l,k) = 0.0;
    
    // Values of rows and columns k and l
    for (uint i: indices) {
        B(i,k) = cosw * aik(i) - sinw * ail(i);
        B(i,l) = cosw * ail(i) + sinw * aik(i);

        R(i,k) = cosw * rik(i) - sinw * ril(i);
        R(i,l) = cosw * ril(i) + sinw * rik(i);
    }

    R(k,k) = cosw * rik(k) - sinw * ril(k);
    R(k,l) = cosw * ril(k) + sinw * rik(k);

    R(l,k) = cosw * rik(l) - sinw * ril(l);
    R(l,l) = cosw * ril(l) + sinw * rik(l);

    B.row(k) = B.col(k).t();
    B.row(l) = B.col(l).t();
    
    *Ev = R;

    return B;
}


double maxoffdiag (arma::mat X, uint *i, uint *j)
{
    // Set diag to zero and find index of extreme value in matrix
    uint ij = arma::abs(arma::trimatu(X) - arma::diagmat(X)).index_max();

    *j = ij / X.n_rows;
    *i = ij - (*j) * X.n_rows;
    
    if (ij == 0) {
        return 0.0;
    } else {
        return std::fabs(X(ij));
    }
}
