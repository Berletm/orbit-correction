#ifndef MODEL_HPP
#define MODEL_HPP

#include "utils.hpp"
#include <vector>

#define G 6.67430e-11

std::vector<Vec3d> compute_accelerations(const std::vector<Object>& objects);

Matrix compute_change_rate(const std::vector<Object>& objects, const Matrix& change_rate);

std::vector<Object> dopri5(std::vector<Object> objects, double dt);

std::vector<Object> integrate(std::vector<Object> objects, double t, double dt);

#endif // MODEL_HPP
