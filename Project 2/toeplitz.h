#ifndef TOEPLITZ_H
#define TOEPLITZ_H

#include <armadillo>

typedef unsigned uint;

arma::mat sym_tridiag (double a, double d, uint N);

class Tridiag
{
    private:
        arma::vec dvar;

        double aconst;
        double dconst;

        uint N;
    public:
        Tridiag (double aconst, double dconst, uint N);

        void set_variable_diags (arma::vec dvar);
        arma::mat get_matrix ();
};

#endif