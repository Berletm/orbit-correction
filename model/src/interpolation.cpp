#include "interpolation.hpp"
#include <stdexcept>
#include <algorithm>

Interpolator::Interpolator(
        const std::vector<Vec3d>& trajectory, 
        const std::vector<double>& times) : times(times)
{
    this->calculate_spline_coef(trajectory);
}

void Interpolator::calculate_spline_coef(const std::vector<Vec3d>& trajectory)
{
    int n = trajectory.size();

    segments.resize(n - 1);
    
    auto solve_for_axis = [&](auto getter) 
    {
        std::vector<double> h(n - 1);
        for (int i = 0; i < n - 1; ++i) h[i] = times[i+1] - times[i];

        std::vector<double> diag(n), upper(n), lower(n), rhs(n);
        std::vector<double> m(n);

        diag[0] = 1.0; upper[0] = 0.0; rhs[0] = 0.0;
        diag[n-1] = 1.0; lower[n-1] = 0.0; rhs[n-1] = 0.0;

        for (int i = 1; i < n - 1; ++i)
        {
            lower[i] = h[i-1];
            diag[i] = 2.0 * (h[i-1] + h[i]);
            upper[i] = h[i];
            
            double f_left = (getter(trajectory[i]) - getter(trajectory[i-1])) / h[i-1];
            double f_right = (getter(trajectory[i+1]) - getter(trajectory[i])) / h[i];
            rhs[i] = 6.0 * (f_right - f_left);
        }

        for (int i = 1; i < n; ++i) 
        {
            double w = lower[i] / diag[i-1];
            diag[i] -= w * upper[i-1];
            rhs[i] -= w * rhs[i-1];
        }

        m[n-1] = rhs[n-1] / diag[n-1];
        for (int i = n - 2; i >= 0; --i) 
        {
            m[i] = (rhs[i] - upper[i] * m[i+1]) / diag[i];
        }
        return m;
    };

    auto mx = solve_for_axis([](const Vec3d& v) { return v.x; });
    auto my = solve_for_axis([](const Vec3d& v) { return v.y; });
    auto mz = solve_for_axis([](const Vec3d& v) { return v.z; });

    for (int i = 0; i < n - 1; ++i) 
    {
        double h = times[i+1] - times[i];
        segments[i].t_start = times[i];
        segments[i].t_end = times[i+1];

        segments[i].ax = trajectory[i].x;
        segments[i].cx = mx[i] / 2.0;
        segments[i].dx = (mx[i+1] - mx[i]) / (6.0 * h);
        segments[i].bx = (trajectory[i+1].x - trajectory[i].x) / h - 
                         (h * (2.0 * mx[i] + mx[i+1])) / 6.0;

        segments[i].ay = trajectory[i].y;
        segments[i].cy = my[i] / 2.0;
        segments[i].dy = (my[i+1] - my[i]) / (6.0 * h);
        segments[i].by = (trajectory[i+1].y - trajectory[i].y) / h - 
                         (h * (2.0 * my[i] + my[i+1])) / 6.0;

        segments[i].az = trajectory[i].z;
        segments[i].cz = mz[i] / 2.0;
        segments[i].dz = (mz[i+1] - mz[i]) / (6.0 * h);
        segments[i].bz = (trajectory[i+1].z - trajectory[i].z) / h - 
                         (h * (2.0 * mz[i] + mz[i+1])) / 6.0;
    }
}

Vec3d Interpolator::interpolate(double time)
{
    auto it = std::upper_bound(times.begin(), times.end(), time);
    int idx;

    if (it == times.begin()) idx = 0;
    else if (it == times.end()) idx = segments.size() - 1;
    else idx = std::distance(times.begin(), it) - 1;

    const auto& seg = segments[idx];
    double dt = time - seg.t_start;
    
    return 
    {
        seg.ax + dt * (seg.bx + dt * (seg.cx + dt * seg.dx)),
        seg.ay + dt * (seg.by + dt * (seg.cy + dt * seg.dy)),
        seg.az + dt * (seg.bz + dt * (seg.cz + dt * seg.dz))
    };
}


