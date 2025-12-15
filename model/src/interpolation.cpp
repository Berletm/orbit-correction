#include "interpolation.hpp"
#include <stdexcept>

Interpolator::Interpolator(
        const std::vector<Vec3d>& trajectory, 
        const std::vector<double>& times) : times(times)
{
    this->calculate_spline_coef(trajectory);
}

void Interpolator::calculate_spline_coef(const std::vector<Vec3d>& trajectory)
{
        int n = trajectory.size();
        
        std::vector<std::vector<double>> derivatives(n, std::vector<double>(3));
        
        for (int i = 0; i < n; i++) 
        {
            if (i == 0) 
            {
                derivatives[i][0] = (trajectory[i + 1].x - trajectory[i].x) / (times[i + 1] - times[i]);
                derivatives[i][1] = (trajectory[i + 1].y - trajectory[i].y) / (times[i + 1] - times[i]);
                derivatives[i][2] = (trajectory[i + 1].z - trajectory[i].z) / (times[i + 1] - times[i]);
            } 
            else if (i == n - 1) 
            {
                derivatives[i][0] = (trajectory[i].x - trajectory[i - 1].x) / (times[i] - times[i - 1]);
                derivatives[i][1] = (trajectory[i].y - trajectory[i - 1].y) / (times[i] - times[i - 1]);
                derivatives[i][2] = (trajectory[i].z - trajectory[i - 1].z) / (times[i] - times[i - 1]);
            } 
            else 
            {
                double dt_left = times[i] - times[i - 1];
                double dt_right = times[i + 1] - times[i];
                
                derivatives[i][0] = ((trajectory[i].x - trajectory[i - 1].x) / dt_left + 
                                    (trajectory[i + 1].x - trajectory[i].x) / dt_right) / 2;
                derivatives[i][1] = ((trajectory[i].y - trajectory[i - 1].y) / dt_left + 
                                    (trajectory[i + 1].y - trajectory[i].y) / dt_right) / 2;
                derivatives[i][2] = ((trajectory[i].z - trajectory[i - 1].z) / dt_left + 
                                    (trajectory[i + 1].z - trajectory[i].z) / dt_right) / 2;
            }
        }
        
        for (int i = 0; i < n - 1; i++) 
        {
            double dt = times[i + 1] - times[i];
            double dt2 = dt * dt;
            double dt3 = dt2 * dt;

            segments[i].ax = trajectory[i].x;
            segments[i].bx = derivatives[i][0];
            segments[i].cx = (3 * (trajectory[i + 1].x - trajectory[i].x) / dt - 
                             2 * derivatives[i][0] - derivatives[i + 1][0]) / dt;
            segments[i].dx = (2 * (trajectory[i].x - trajectory[i + 1].x) / dt + 
                             derivatives[i][0] + derivatives[i + 1][0]) / dt2;
            
            segments[i].ay = trajectory[i].y;
            segments[i].by = derivatives[i][1];
            segments[i].cy = (3 * (trajectory[i + 1].y - trajectory[i].y) / dt - 
                             2 * derivatives[i][1] - derivatives[i + 1][1]) / dt;
            segments[i].dy = (2 * (trajectory[i].y - trajectory[i + 1].y) / dt + 
                             derivatives[i][1] + derivatives[i + 1][1]) / dt2;
            
            segments[i].az = trajectory[i].z;
            segments[i].bz = derivatives[i][2];
            segments[i].cz = (3 * (trajectory[i + 1].z - trajectory[i].z) / dt - 
                             2 * derivatives[i][2] - derivatives[i + 1][2]) / dt;
            segments[i].dz = (2 * (trajectory[i].z - trajectory[i + 1].z) / dt + 
                             derivatives[i][2] + derivatives[i + 1][2]) / dt2;
        }
    }

Vec3d Interpolator::interpolate(double time)
{   
    if (time <= times.front()) 
    {
        Vec3d res = {segments[0].ax, segments[0].ay, segments[0].az};
        return res;
    }

    if (time >= times.back()) 
    {
        int last = segments.size() - 1;
        double dt = time - segments[last].t_start;

        Vec3d res = 
        {
            segments[last].ax + segments[last].bx * dt + segments[last].cx * dt * dt + segments[last].dx * dt * dt * dt,
            segments[last].ay + segments[last].by * dt + segments[last].cy * dt * dt + segments[last].dy * dt * dt * dt,
            segments[last].az + segments[last].bz * dt + segments[last].cz * dt * dt + segments[last].dz * dt * dt * dt
        };

        return res;
    }

    for (const auto& seg : segments) 
    {
        if (time >= seg.t_start && time <= seg.t_end) 
        {
            double dt = time - seg.t_start;
            double dt2 = dt * dt;
            double dt3 = dt2 * dt;
            
            Vec3d res =
            {
                seg.ax + seg.bx * dt + seg.cx * dt2 + seg.dx * dt3,
                seg.ay + seg.by * dt + seg.cy * dt2 + seg.dy * dt3,
                seg.az + seg.bz * dt + seg.cz * dt2 + seg.dz * dt3
            };

            return res;
            
        }
    }

    throw std::runtime_error("Time out of range");
}


