#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <string>
#include <time.h>

int main(int argc, char **argv)
{
    double relative_max_error (arma::vec u, arma::vec v);
    arma::vec f (arma::vec x); 
    arma::vec u_exact (arma::vec x);
    arma::vec set_and_fill (unsigned int n, double filler); 
    arma::vec poisson_solver (arma::vec f_vec, unsigned int n);
    arma::vec tridiag_solver (arma::vec f_vec, arma::vec a, arma::vec b,
                              arma::vec c, unsigned int n);

    arma::vec n_vec; 
    n_vec << 10 << 100 << 1000 << 10000 << 100000 << 
             1000000 << 10000000;
    
    arma::vec re_max(n_vec.n_elem);
    arma::vec h_vec(n_vec.n_elem);
    
    arma::vec f_vec;
    arma::vec v;
    arma::vec u;
    arma::vec x;

    arma::vec a;
    arma::vec b;
    arma::vec c;

    // Filenames for u and v text files    
    std::string nstr;
    std::string ustr;
    std::string vstr;

    // File for comparing poisson and general tridiagonal method
    std::ofstream method_compare;

    clock_t start;
    clock_t end;

    double general_time;
    double poisson_time;
    double hh;
    double h;

    unsigned int n;
    unsigned int k = 4; // For getting only wanted values of n_vec in the loop 
    unsigned int time_runs = 10;
    
    for (int i = 0; i < n_vec.n_elem - k; i++) {
        n = n_vec[i];
        h = 1. / ((double) n + 1);
        hh = h * h;
        h_vec[i] = h;

        // Fills vectors with a value
        a = set_and_fill(n, -1.0);
        b = set_and_fill(n, 2.0);
        c = a;
        
        x = arma::linspace<arma::vec>(h, 1.0 - h, n);
        u = u_exact(x);
        f_vec = hh * f(x);
        
        for (int j = 0; j < time_run; j++) {
            start = clock();
            v = tridiag_solver(f_vec, a, b, c, n);
            end = clock();
            general_time += (double) (end - start) / CLOCKS_PER_SEC;
        }
        general_time /= (double) time_run;

        // Saving re_max for later analysis
        re_max[i] = relative_max_error(u, v);

        // Filenames for saved vectors
        nstr = std::to_string(n) + ".txt";
        ustr = "results/u" + nstr;
        vstr = "results/v" + nstr;

        u.save(ustr, arma::raw_ascii);
        v.save(vstr, arma::raw_ascii);

        // Comaparing timeing for 
    }

    method_compare.open("results/method_compare.txt");
    
    for (int i = n_vec.n_elem - k; i < n_vec.n_elem; i++) {
        method_compare << "------------------------------" << std::endl;

        n = n_vec[i];
        h = 1.0 / (n + 1.0);
        hh = h * h;
        h_vec[i] = h;
        
        x = arma::linspace<arma::vec>(h, 1.0 - h, n);
        u = u_exact(x);
        f_vec = hh * f(x);

        a = set_and_fill(n, -1.0);
        b = set_and_fill(n, 2.0);
        c = a;
        
        method_compare << "N = " << n << std::endl;

        // Timing the general tridiagonal solver
        general_time = 0;
        for (int j = 0; j < time_runs; j++) {
            start = clock();
            v = tridiag_solver(f_vec, a, b, c, n);
            end = clock();
            general_time += (double) (end - start) / CLOCKS_PER_SEC;
        }
        general_time /= (double) time_runs;

        method_compare << "General solver = " << general_time <<
                          " seconds." << std::endl;

        // Saving maximum relative error for this n
        re_max[i] = relative_max_error(u, v);
        
        // Reset v to make sure performance not impacted
        v.set_size(n);
       
        // Timing the poisson special solver
        poisson_time = 0;
        for (int j = 0; j < time_runs; j++) {
            start = clock();
            v = poisson_solver(f_vec, n);
            end = clock();
            poisson_time += (double) (end - start) / CLOCKS_PER_SEC;
        }
        poisson_time /= (double) time_runs;

        method_compare << "Poisson solver = " << poisson_time <<
                          " seconds." << std::endl;

        method_compare << "Time of poisson solver / " << 
                          "Time of general solver = " <<
                          poisson_time / general_time << std::endl;
    }
    
    method_compare.close();

    re_max = arma::log10(re_max);
    re_max.save("results/re_max.txt", arma::raw_ascii);
    h_vec.save("results/h.txt", arma::raw_ascii);



    return 0;
}


arma::vec tridiag_solver(arma::vec f_vec, arma::vec a, arma::vec b,
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


arma::vec f(arma::vec x)
{
    return 100.0 * arma::exp(-10.0 * x);
}


arma::vec u_exact(arma::vec x)
{
    return 1.0 - (1.0 - std::exp(-10.0)) * x - arma::exp(-10 * x);
}


double relative_max_error(arma::vec u, arma::vec v) 
{
    arma::vec re = arma::abs((u - v) / u);
    return re.max();
}


arma::vec set_and_fill(unsigned int n, double filler)
{
    return arma::ones<arma::vec>(n) * filler;
}
