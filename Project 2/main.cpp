#include <iostream>

#include "jacobi.h"

int main(int argc, char **argv)
{
    arma::arma_rng::set_seed_random();

    arma::mat A(10, 10);
    A.randu();
    A = A.t() * A;

    arma::vec eigval = arma::eig_sym(A);

    A.print();
    A = jacobi(A, 1.0E-16, 10000);
    std::cout << "-----------------------------------" << std::endl;
    A.diag().print();

    eigval.print();
    return 0;
}