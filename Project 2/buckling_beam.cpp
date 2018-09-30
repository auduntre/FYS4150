#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "jacobi.h"
#include "toeplitz.h"


int main(int argc, char **argv)
{
    void set_values(double *h, double *a, double *d, double *deps, 
                    double eps, uint N);
    bool obtained_analy(arma::vec val, arma::vec lambda, 
                        std::string method_name, double deps); 
    arma::vec analy_eigval(double a, double d, uint N);
    

    arma::vec eigval;
    arma::vec jacobival;
    arma::vec lambda;
    
    arma::mat A;
    arma::mat Ev;

    double eps = 1.0E-12;
    double deps;
    double h;
    double a;
    double d;

    uint Ns[3]= {10, 50, 100};

    for (uint N: Ns) {
        std::cout << "--------------------" << std::endl;
        std::cout << "N = " << N << std::endl;

        set_values(&h, &a, &d, &deps, eps, N);
        Tridiag td(a, d, N-1);

        A = td.get_matrix(); //sym_tridiag (a, d, N-1);
        Ev = arma::eye(A.n_rows, A.n_cols);
        
        lambda = analy_eigval(a, d, N);
        eigval = arma::eig_sym(A);
        jacobival = arma::sort(jacobi(A, &Ev, eps));

        if (!obtained_analy(eigval, lambda, "Amadillo", deps) ||
            !obtained_analy(jacobival, lambda, "Jacobi", deps)) {
            return 1;
        }
    }
    return 0;
}


void set_values(double *h, double *a, double *d, double *deps, 
                double eps, uint N) 
{
    double delta = 1.0 / (double) N ;

    *deps = eps * N * 10;
    *h = delta;
    *a = -1.0 / (delta * delta);
    *d = 2.0 / (delta * delta);
}


bool obtained_analy(arma::vec val, arma::vec lambda, std::string method_name,
                    double deps) {
    if (arma::approx_equal(val, lambda, "absdiff", deps)) {
        std::cout << method_name << " obtained approximately "
                  << "analytical eigenvalues with eps = " 
                  << deps << std::endl;
        return true;
    } else {
        std::cout << (val - lambda).max() << std::endl;
        return false;
    }
}


arma::vec analy_eigval(double a, double d, uint N)
{
    arma::vec js = arma::linspace<arma::vec>(1, N-1, N-1);
    arma::vec lambda(N-1);

    return d + 2 * a * arma::cos(js * arma::datum::pi/ N);
}