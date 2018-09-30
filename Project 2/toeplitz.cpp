#include "toeplitz.h"
#include <iostream>

Tridiag::Tridiag (double aconst, double dconst, uint N)
{
    this->N = N;
    
    this->dvar = arma::zeros<arma::vec>(this->N);    
    this->aconst = aconst;
    this->dconst = dconst;
}


void Tridiag::set_variable_diags (arma::vec dvar)
{
    for (uint i = 0; i < this->dvar.n_elem; i++) {
        this->dvar(i) = dvar(i); 
    }
}


arma::mat Tridiag::get_matrix ()
{
    arma::vec diagc = arma::zeros<arma::vec>(this->N);
    
    diagc(0) = this->dconst ;
    diagc(1) = this->aconst ;

    return arma::toeplitz(diagc) + arma::diagmat(this->dvar);
}