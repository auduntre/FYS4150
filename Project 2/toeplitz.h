#ifndef TOEPLITZ_H
#define TOEPLITZ_H

#include <armadillo>

typedef unsigned uint;


/*
 * Class for setting up an tridiagonal Toeplitz' matrix (or one
 * with a variable diagonal) with the same value over and under 
 * the diagonal.
 */
class Tridiag
{
    private:
        arma::vec dvar;

        double aconst;
        double dconst;

        uint N;
    public:
        /*
         * Setting up tridiagonal Toeplitz' matrix in Armadillo
         * with the constant elements on the diagonal and under/over it.
         *
         * @param aconst: Value of the elements over/under the diagonal.
         * @param dconst: The constant value of the diagonal
         * @param N: NxN is the size of the Armadillo matrix.
         */
        Tridiag (double aconst, double dconst, uint N);

        /*
         * Setting up variable elements on the diagonal
         *
         * @param dvar: Armadillo vector with the variable elements.
         */
        void set_variable_diags (arma::vec dvar);

        /*
         * Get method for the tridiagonal Toeplitz matrix.
         *
         * @return The tridiagonal Toeplitz, Armadillo matrix.
         */
        arma::mat get_matrix ();
};

#endif
