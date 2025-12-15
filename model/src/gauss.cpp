#include "gauss.hpp"

Solver::Solver(const Matrix& A, const Matrix& b) : A(A), b(b) {}

Matrix Solver::solve()
{
    forward_elimination();    

    return backward_substitution();
}

void Solver::swap_row(int i, int j)
{
    for (int k = 0; k < A.cols; ++k)
    {
        double temp = A.mat[i][k];
        A.mat[i][k] = A.mat[j][k];
        A.mat[j][k] = temp;
    }
}

void Solver::forward_elimination()
{
    for (int k = 0; k < A.rows; ++k)
    {
        int cur_max_idx = k;
        double cur_max  = A.mat[cur_max_idx][k];
        
        for (int i = k + 1; i < A.cols; ++i)
        {
            if (std::abs(A.mat[i][k]) > cur_max) cur_max = A.mat[i][k], cur_max_idx = i;
        }

        if (cur_max_idx != k) swap_row(k, cur_max_idx);

        for (int i = k + 1; i < A.rows; ++i)
        {
            double factor = A.mat[i][k] / A.mat[k][k];
            
            for (int j = k + 1; j < A.rows; ++j)
            {
                A.mat[i][j] -= A.mat[k][j] * factor;
                b.mat[j][0] -=  b.mat[k][0] * factor;
            }
            A.mat[i][k] = 0;
        }
    }
}

Matrix Solver::backward_substitution()
{
    Matrix res(b.rows, b.cols);

    for (int i = A.rows - 1; i >= 0; --i)
    {
        res.mat[i][0] = b.mat[i][0];

        for (int j = i + 1; j < A.cols; ++j)
        {
            res.mat[i][0] -= A.mat[i][j] * res.mat[j][0];
        }

        res.mat[i][0] = res.mat[i][0] / A.mat[i][i];
    }

    return res;
}