#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <string>


int main(int argc, char **argv)
{
    arma::vec poisson_solver (arma::vec f_vec, unsigned int n);
    arma::vec solver (arma::vec f_vec, arma::vec a, arma::vec b,
                      arma::vec c, unsigned int n);
    double u_exact (double x);
    double f (double x); 

    arma::vec n_vec;
    n_vec << 10 << 100 << 1000 << 1000000 << 10000000;
    
    std::string ext = ".txt";
    std::string nstr;
    std::string ustr;
    std::string vstr;

    double hh;
    double h;
    int n;

    for (int i = 0; i < n_vec.n_elem - 2; i++) {

        n = n_vec[i];
        h = 1. / ((double) n + 1);
        hh = h * h;

        arma::vec x = arma::linspace<arma::vec> (0.0, 1.0, n+2);
        
        arma::vec f_vec(n);
        arma::vec a(n);
        arma::vec b(n);
        arma::vec c(n);
        
        a.fill(-1.0);
        b.fill(2.0);
        c.fill(-1.0);
        
        arma::vec v(n);
        arma::vec u(n);
        arma::vec re(n);

        for (int j = 0; j < n; j++) {
            u[j] = u_exact(x[j+1]);
            f_vec[j] = f(x[j+1]) * hh; 
        }

        //v = solver(f_vec, a, b, c, n);
        v = poisson_solver(f_vec, n);

        for (int j = 0; j < n; j++) {
            re[j] = (u[j] - v[j]) / u[j];
        }

        nstr = std::to_string(n) + ext;
        ustr = "u" + nstr;
        vstr = "v" + nstr;

        u.save(ustr, arma::raw_ascii);
        v.save(vstr, arma::raw_ascii);
    }

    return 0;
}


arma::vec solver(arma::vec f_vec, arma::vec a, arma::vec b,
                 arma::vec c, unsigned int n)
{
    arma::vec v(n);

    //Forward subst
    for (int i = 1; i < n; i++) {
        b[i] = b[i] - (a[i-1] * c[i-1]) / b[i-1];
        f_vec[i] = f_vec[i] - (a[i-1] * f_vec[i-1]) / b[i-1];
    }

    //Backward subst
    v[n-1] = f_vec[n-1] / b[n-1];
    for (int i = n-2; i >= 0; i--) {
        v[i] = (f_vec[i] - c[i] * v[i+1]) / b[i];
    }

    return v;
}


arma::vec poisson_solver(arma::vec f_vec, unsigned int n)
{
    arma::vec v(n);
    arma::vec b(n);

    b[0] = 2.0; 
    for (int i = 1; i < n; i++) {
        b[i] = 2.0 - 1 / b[i-1];
        f_vec[i] = f_vec[i] + (f_vec[i-1] / b[i-1]);
    }

    v[n-1] = f_vec[n-1] / b[n-1]; 
    for (int i = n-2; i >= 0; i--) {
        v[i] = (f_vec[i] + v[i+1]) / b[i];
    }
    
    return v;
}


double f(double x)
{
    return 100.0 * std::exp(-10.0 * x);

}


double u_exact(double x)
{
    return 1.0 - (1.0 - std::exp(-10.0)) * x - std::exp(-10 * x);
}