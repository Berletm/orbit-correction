#ifndef GAUSS_ELIMINATION_HPP
#define GAUSS_ELIMINATION_HPP

#include "utils.hpp"

class Solver
{
private:
    Matrix A, b;

    void forward_elimination();
    
    Matrix backward_substitution();

    void swap_row(int i, int j);
public:
    Solver(const Matrix& A, const Matrix& b);

    Matrix solve();
};


#endif // GAUSS_ELIMINATION_HPP
