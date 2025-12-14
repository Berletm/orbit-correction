#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "utils.hpp"
#include <vector>

class Interpolator
{
private:
    struct SplineSegment
    {
        double ax, bx, dx, cx;
        double ay, by, dy, cy;
        double az, bz, dz, cz;
        double t_start, t_end;
    };

    std::vector<SplineSegment> segments;
    std::vector<double> times;
    
    void calculate_spline_coef(const std::vector<Vec3d>& trajectory);
public:
    Interpolator(
        const std::vector<Vec3d>& trajectory, 
        const std::vector<double>& times);

    Vec3d interpolate(double time);
};


#endif // INTERPOLATION_HPP
