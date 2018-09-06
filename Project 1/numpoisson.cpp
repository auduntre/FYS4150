#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <string>
#include <time.h>

int main(int argc, char **argv)
{
    arma::vec set_and_fill (unsigned int n, double filler);
    arma::vec poisson_solver (arma::vec f_vec, unsigned int n);
    arma::vec solver (arma::vec f_vec, arma::vec a, arma::vec b,
                      arma::vec c, unsigned int n);
    double u_exact (double x);
    double f (double x); 

    arma::vec n_vec; 
    arma::vec re_max(n_vec.n_elem);
    arma::vec re;
    arma::vec f_vec;
    arma::vec v;
    arma::vec u;
    arma::vec x;

    arma::vec a;
    arma::vec b;
    arma::vec c;
    
    std::string nstr;
    std::string ustr;
    std::string vstr;

    std::ofstream method_compare;

    clock_t start;
    clock_t end;
    
    double method_time;
    double hh;
    double h;
    int n;
    
    n_vec << 10 << 100 << 1000 << 1000000 << 10000000;
    
    for (int i = 0; i < n_vec.n_elem - 2; i++) {
        n = n_vec[i];
        h = 1. / ((double) n + 1);
        hh = h * h;

        x = arma::linspace<arma::vec>(0.0, 1.0, n+2);
        
        a = set_and_fill(n, -1.0);
        b = set_and_fill(n, 2.0);
        c = a;
        
        f_vec.set_size(n);
        u.set_size(n);
        re.set_size(n);

        for (int j = 0; j < n; j++) {
            u[j] = u_exact(x[j+1]);
            f_vec[j] = f(x[j+1]) * hh; 
        }

        v = solver(f_vec, a, b, c, n);

        // Saving re_max for later analysis
        for (int j = 0; j < n; j++) {
            re[j] = std::fabs((u[j] - v[j]) / u[j]);
        }
        re_max[i] = re.max();

        nstr = std::to_string(n) + ".txt";
        ustr = "u" + nstr;
        vstr = "v" + nstr;

        u.save(ustr, arma::raw_ascii);
        v.save(vstr, arma::raw_ascii);
    }

    method_compare.open("method_compare.txt");
    
    for (int i = n_vec.n_elem - 3; i < n_vec.n_elem; i++) {
        n = n_vec[i];
        x = arma::linspace<arma::vec>(0.0, 1.0, n+2);

        v.set_size(n);
        f_vec.set_size(n);
        a = set_and_fill(n, -1.0);
        b = set_and_fill(n, 2.0);
        c = set_and_fill(n, -1.0);
        
        h = 1.0 / (n + 1.0);
        hh = h * h;

        method_compare << "N = " << n << std::endl;

        for (int  j = 0; j < n; j++) {
            f_vec[j] = f(x[j+1]) * hh;
        }

        start = clock();
        v = solver(f_vec, a, b, c, n);
        end = clock();

        method_time = (end - start) / CLOCKS_PER_SEC;
        method_compare << "General solver = " << method_time <<
                          " seconds." << std::endl;

        v.set_size(n);
        
        start = clock();
        v = poisson_solver(f_vec, n);
        end = clock();

        method_time = 1000 * (end - start) / CLOCKS_PER_SEC;
        method_compare << "Poisson solver = " << method_time <<
                          " seconds." << std::endl;
       
    }

    method_compare.close();
    
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

    //Forward subst
    b[0] = 2.0; 
    for (int i = 1; i < n; i++) {
        b[i] = 2.0 - 1 / b[i-1];
        f_vec[i] = f_vec[i] + (f_vec[i-1] / b[i-1]);
    }

    //Backward subst
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


arma::vec set_and_fill(unsigned int n, double filler)
{
    arma::vec y(n);
    return y.fill(filler);
}
