#ifndef SOLVE_TRI_H
#define SOLVE_TRI_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace spline_interpolation
{
using TriDiagonal = Eigen::MatrixXd;
using Diagonal = Eigen::VectorXd;

/* The following function implements the Thomas algorithm for solving tridiagonal linear systems of equations
 * in the form Ax=B, where A the tridiagonal matrix and B the values matrix.
 *
 * The function checks the sizes of the given matrices and will throw an exception if they don't much with the
 * expected ones:
 * - Tridiagonal matrix A -> N x N
 * - Values matrix B -> N x M
 *
 * Returns x where x is a N x M matrix that satisfies the linear system.
 *
 * The algorithm does not modify the objects given in the arguments of the function.
 * For more info: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 */
Eigen::MatrixXd solveTridiagonal(const TriDiagonal& tridiagonal_matrix, Eigen::MatrixXd values_matrix) {
    // Check if tridiagonal matrix is square and if the values matrix size matches as well
    const int matrix_rows{ static_cast<int>(tridiagonal_matrix.rows()) };

    if (matrix_rows != tridiagonal_matrix.cols())
    {
        // The tridiagonal matrix is not actually square and an error is thrown
        throw std::runtime_error("SolveTridiagonal() -> Given matrix is not square");
    }
    if (matrix_rows != values_matrix.rows())
    {
        // The size of the values matrix does not match with the tridiagonal matrix and an error is thrown
        throw std::runtime_error("SolveTridiagonal() -> Given values matrix does not match tridiagonal");
    }

    // Create the result matrix based on the dimensions of the given matrices
    Eigen::MatrixXd return_matrix{ matrix_rows, values_matrix.cols() };

    // Get the three diagonals
    Diagonal main_diagonal{ tridiagonal_matrix.diagonal(0) };
    Diagonal super_diagonal{ tridiagonal_matrix.diagonal(1) };
    Diagonal sub_diagonal{ tridiagonal_matrix.diagonal(-1) };

    // Do Forward Sweep Steps
    for (int i{ 1 }; i < matrix_rows; i++)
    {
        double temp{ sub_diagonal(i - 1) / main_diagonal(i - 1) };
        main_diagonal(i) = main_diagonal(i) - temp * super_diagonal(i - 1);
        values_matrix.row(i) = values_matrix.row(i) - temp * values_matrix.row(i - 1);
    }
    // Do Back Substituion
    return_matrix.row(matrix_rows - 1) = values_matrix.row(matrix_rows - 1) / main_diagonal(matrix_rows - 1);
    for (int i{ matrix_rows - 2 }; i >= 0; i--)
    {
        return_matrix.row(i) = (values_matrix.row(i) - super_diagonal(i) * return_matrix.row(i + 1)) / main_diagonal(i);
    }

    return return_matrix;
}
} // namespace spline_interpolation
#endif