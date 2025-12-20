#ifndef CORRECTION_HPP
#define CORRECTION_HPP

#include <vector>
#include "utils.hpp"

std::vector<Celestial> calculate_residuals(const std::vector<Celestial>& observed, const std::vector<Celestial>& computed);

std::vector<Celestial> cart2celestial(const std::vector<Vec3d>& coords);

Matrix calculate_state_change(const Vec3d& computed);

Matrix calculate_jacobian(const Matrix& state_change, const Matrix& change_rate);

Matrix solve(const Matrix& A, const Matrix& b);

Matrix stack_matrix(const std::vector<Matrix>& matrices);

Matrix stack_vector(const std::vector<Celestial>& vecs);

SystemState interpolate(
    const std::vector<std::vector<Object>>& object_trajectories, 
    const std::vector<SystemState>& trajectory, 
    double t);

void correction(
    std::vector<Object>& initial_state, 
    const std::vector<Celestial>& observed, 
    const std::vector<Vec3d>& obs_position, 
    const std::vector<double>& obs_time,
    const Matrix& weights
);

void read_observed_data(
    std::vector<double>& time, 
    std::vector<Vec3d>& observatories, 
    std::vector<Celestial>& observed,
    Matrix& weights    
);


#endif //CORRECTION_HPP
