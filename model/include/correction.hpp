#ifndef CORRECTION_HPP
#define CORRECTION_HPP

#include <vector>
#include "utils.hpp"

std::vector<Celestial> calculate_residuals(const std::vector<Celestial>& observed, const std::vector<Celestial>& computed);

std::vector<Celestial> cart2celestial(const std::vector<Vec3d>& coords);

Matrix compute_change_rate();


#endif //CORRECTION_HPP
