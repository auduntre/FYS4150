#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <cmath>

typedef unsigned uint;

/*  
 *  Jacobi's rotational method for diagonalizing a symmetric matrix
 *  untill all the off-diagonal elements are less the a given tolerance, 
 *  then returns the now eigenvalues on the diagonal as an vector.
 *
 *  @param A: symmetric Armadillo matrix to be diagonlized.
 *  @param *Ev: Pointer to an identity Armadillo matrix of 
 *  same dimensions as A. Will be transformed into an matrix 
 *  with the eigenvectors of A as column.
 *  @param eps: Tolerance for the off-diagonal elements.
 *  @return An Armadillo vector with the eigenvalues of A.
 */ 
arma::vec jacobi (arma::mat A, arma::mat *Ev, double eps);


/*
 * Jacobi rotation of a matrix on the rows and columns given.
 *
 * @param B: Symmetric Armadillo matrix to be rotated.
 * @param *Ev: Pointer to Armadillo matrix for eigenvectors of B. 
 * @param k: Row and column to be rotated, k < l.
 * @param l: Row and column to be rotated, l > k.
 * @return the rotated matrix B.
 */
arma::mat rotate (arma::mat B, arma::mat *Ev, uint k, uint l);


/*
 * Finds the maximum off-diagonal element in absolute value of
 * a symmetric matrix.
 *
 * @param X. Symmetric Amadillo matrix
 * @param *i: Pointer to the maximum off-diagonal element's row.
 * @param *j: Pointer to the maximum off-diagonal element's column.
 * @return the absolute value of the off-diagonal element.
 */
double maxoffdiag (arma::mat X, uint *i, uint *j);

#endif
