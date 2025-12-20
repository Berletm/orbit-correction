#ifndef MODEL_HPP
#define MODEL_HPP

#include "utils.hpp"
#include <vector>

#define G 6.67430e-11

std::vector<Vec3d> compute_accelerations(const std::vector<Object>& objects);

Matrix compute_gravitational_gradient(const std::vector<Object>& objects);

Matrix compute_change_rate(const std::vector<Object>& objects, const Matrix& change_rate);

std::vector<Object> dopri5(std::vector<Object> objects, SystemState& state, double dt);

void integrate(
    std::vector<Object> objects,
    std::vector<SystemState>& states,
    std::vector<std::vector<Object>>& objects_trajectories,
    double t, double dt);

#endif // MODEL_HPP
