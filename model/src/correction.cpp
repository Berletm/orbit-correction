#include "correction.hpp"
#include <cmath>

std::vector<Celestial> calculate_residuals(const std::vector<Celestial>& observed, const std::vector<Celestial>& computed)
{
    std::vector<Celestial> res;
    for (int i = 0; i < observed.size(); ++i)
    {
        Celestial delta = observed[i] - computed[i];
        res.push_back(delta);
    }

    return res;
}

std::vector<Celestial> cart2celestial(const std::vector<Vec3d>& coords)
{
    std::vector<Celestial> res;

    for (int i = 0; i < coords.size(); ++i)
    {
        Vec3d current_coords = coords[i];
        
        double ra  = std::atan2(current_coords.y, current_coords.x);
        double dec = std::atan2(current_coords.z, std::sqrt(current_coords.y * current_coords.y + current_coords.x * current_coords.x));

        Celestial current_angles(ra, dec);
        res.push_back(current_angles);
    }
    
    return res;
}
