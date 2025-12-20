#include "correction.hpp"
#include "model.hpp"
#include "interpolation.hpp"
#include "gauss.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

double angle2rad(double angle)
{
    return angle / 180 * M_PI;
}

double calculate_error(const std::vector<Celestial>& residuals)
{
    double ans = 0.0;

    for (const auto& res: residuals)
    {
        ans += (std::pow(res.ra, 2) + std::pow(res.dec, 2));
    }

    return ans;
}

std::vector<Celestial> calculate_residuals(const std::vector<Celestial>& observed, const std::vector<Celestial>& computed)
{
    std::vector<Celestial> res;
    for (int i = 0; i < observed.size(); ++i)
    {
        Celestial delta = observed[i] - computed[i];

        while (delta.ra >  M_PI) delta.ra -= 2 * M_PI;
        while (delta.ra < -M_PI) delta.ra += 2 * M_PI;

        while (delta.dec >  M_PI) delta.dec -= 2 * M_PI;
        while (delta.dec < -M_PI) delta.dec += 2 * M_PI;

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

Matrix calculate_state_change(const Vec3d& computed)
{
    Matrix res(2, 3);

    double dec = std::atan2(computed.z, std::sqrt(computed.y * computed.y + computed.x * computed.x));

    double general = pow(computed.x, 2) + pow(computed.y, 2);
    double dist_sq = computed * computed;
    double cos_dec = cos(dec);
    res.mat[0][0] = - (computed.y / general) * cos_dec;
    res.mat[0][1] = (computed.x / general) * cos_dec;
    res.mat[0][2] = 0;

    res.mat[1][0] = - computed.x * computed.z / (sqrt(general) * dist_sq);
    res.mat[1][1] = - computed.y * computed.z / (sqrt(general) * dist_sq);
    res.mat[1][2] = sqrt(general) / dist_sq;

    return res;
}

Matrix calculate_jacobian(const Matrix& state_change, const Matrix& change_rate)
{
    Matrix j = state_change * change_rate * -1;

    return j;
}

void read_observed_data(std::vector<double>& time, std::vector<Vec3d>& observatories, std::vector<Celestial>& observed)
{
    std::ifstream file("../data/output.txt");
    
    if (!file.is_open()) return;

    std::string line;

    while (std::getline(file, line))
    {
        std::istringstream ss(line);

        double t = 0;
        Vec3d observatory(0, 0, 0);
        Celestial obs(0, 0);

        ss >> t >> obs.ra >> obs.dec >> observatory.x >> observatory.y >> observatory.z;

        obs.ra = angle2rad(obs.ra);
        obs.dec = angle2rad(obs.dec);

        time.push_back(t);
        observatories.push_back(observatory);
        observed.push_back(obs);
    }
    
    file.close();
}

Matrix stack_matrix(const std::vector<Matrix>& matrices)
{
    Matrix res(matrices[0].rows * matrices.size(), matrices[0].cols);

    for (int idx = 0; idx < matrices.size(); ++idx)
    {
        Matrix mat = matrices[idx];
        int shift = idx * mat.rows;

        for (int i = 0; i < mat.rows; ++i)
        {
            for (int j = 0; j < mat.cols; ++j)
            {
                res.mat[shift + i][j] = mat.mat[i][j];
            }
        }
    }

    return res;
}

Matrix stack_vector(const std::vector<Celestial>& vecs)
{
    Matrix res(2 * vecs.size(), 1);

    for (int i = 0; i < vecs.size(); ++i)
    {   
        int shift = i * 2;
        Celestial cur_vec = vecs[i];

        res.mat[0 + shift][0] = cur_vec.ra * cos(cur_vec.dec); 
        res.mat[1 + shift][0] = cur_vec.dec;
    }

    return res;
}

Matrix solve(const Matrix& A, const Matrix& b)
{
    Solver solver(A, b);

    return solver.solve();
}

Matrix interpolate_change_rate(const std::vector<Matrix>& trajectory, const std::vector<double>& times, double target_time)
{
    std::vector<Vec3d> row_1, row_2, row_3;

    for(const auto& mat: trajectory)
    {
        row_1.push_back({mat.mat[0][0], mat.mat[0][1], mat.mat[0][2]});
        row_2.push_back({mat.mat[1][0], mat.mat[1][1], mat.mat[1][2]});
        row_3.push_back({mat.mat[2][0], mat.mat[2][1], mat.mat[2][2]});
    }

    Interpolator row_1_interpolator(row_1, times);
    Interpolator row_2_interpolator(row_2, times);
    Interpolator row_3_interpolator(row_3, times);

    Vec3d row_1_target = row_1_interpolator.interpolate(target_time);
    Vec3d row_2_target = row_2_interpolator.interpolate(target_time);
    Vec3d row_3_target = row_3_interpolator.interpolate(target_time);

    Matrix res(3, 3);

    res.mat[0][0] = row_1_target.x;
    res.mat[0][1] = row_1_target.y;
    res.mat[0][2] = row_1_target.z;

    res.mat[1][0] = row_2_target.x;
    res.mat[1][1] = row_2_target.y;
    res.mat[1][2] = row_2_target.z;

    res.mat[2][0] = row_3_target.x;
    res.mat[2][1] = row_3_target.y;
    res.mat[2][2] = row_3_target.z;

    return res;
}

SystemState interpolate(const std::vector<std::vector<Object>>& object_trajectories, const std::vector<SystemState>& trajectory, double t)
{
    std::vector<Vec3d> oumuamua_trajectory;
    std::vector<Vec3d> earth_trajectory;
    std::vector<Vec3d> sun_trajectory;
    std::vector<Vec3d> jupiter_trajectory;
    std::vector<Matrix> rate_change;
    std::vector<double> times;

    for (const auto& state: trajectory)
    {
        oumuamua_trajectory.push_back(state.positions[0]);
        sun_trajectory.push_back(state.positions[1]);
        jupiter_trajectory.push_back(state.positions[2]);
        earth_trajectory.push_back(state.positions[3]);
        times.push_back(state.time);
        rate_change.push_back(state.change_rate);
    }

    Interpolator oumuamua_interpolator(oumuamua_trajectory, times);
    Interpolator earth_interpolator(earth_trajectory, times);
    Interpolator sun_interpolator(sun_trajectory, times);
    Interpolator jupiter_interpolator(jupiter_trajectory, times);

    SystemState res;

    res.positions.push_back(oumuamua_interpolator.interpolate(t));
    res.positions.push_back(sun_interpolator.interpolate(t));
    res.positions.push_back(jupiter_interpolator.interpolate(t));
    res.positions.push_back(earth_interpolator.interpolate(t));

    res.change_rate = interpolate_change_rate(rate_change, times, t);

    return res;
}

Vec3d mat2vec(const Matrix& mat)
{
    Vec3d res;

    res.x = mat.mat[0][0];
    res.y = mat.mat[1][0];
    res.z = mat.mat[2][0];

    return res;
}

void write_residuals(const std::vector<Celestial>& residuals, std::ofstream& file)
{
    file << std::setprecision(15);

    for (const auto& r: residuals)
    {
        file << r.ra << " " << r.dec << "\n";
    }
    
    file << "squared error: " << calculate_error(residuals) << "\n";
    file << "----------------------\n";
}

void correction(std::vector<Object>& initial_state, const std::vector<Celestial>& observed, const std::vector<Vec3d>& obs_position, const std::vector<double>& obs_time)
{   
    std::ofstream file("residuals.txt");
    Vec3d params = initial_state[0].position;

    double t_end = *std::max_element(obs_time.begin(), obs_time.end());

    for (int i = 0; i < 5; ++i)
    {   
        std::vector<Object> initial_state_copy = initial_state;
        initial_state_copy[0].position = params;

        std::vector<SystemState> states;
        std::vector<std::vector<Object>> objects_trajectories;
        integrate(initial_state_copy, states, objects_trajectories, t_end, 60*60*12); // step = 12 hours

        std::vector<Vec3d> computed_positions;
        std::vector<Matrix> computed_jacobians;

        for (int i = 0; i < obs_time.size(); ++i)
        {
            double t = obs_time[i];
            SystemState state = interpolate(objects_trajectories, states, t);
            Vec3d earth_pos = state.positions[3]; 
            Vec3d oumuamua_pos = state.positions[0];

            Vec3d rho_vec = oumuamua_pos - (earth_pos + obs_position[i]);
            computed_positions.push_back(rho_vec);

            Matrix state_change = calculate_state_change(rho_vec);
            Matrix j = calculate_jacobian(state_change, state.change_rate);
            computed_jacobians.push_back(j);
        }

        std::vector<Celestial> computed_angles = cart2celestial(computed_positions);

        std::vector<Celestial> r = calculate_residuals(observed, computed_angles);
        write_residuals(r, file);

        Matrix J = stack_matrix(computed_jacobians); // 2N x 3 (rows x cols)
        Matrix R = stack_vector(r); // 2N x 1 (rows x cols)

        Matrix Jt = J.transposed(); // 3 x 2N
        Matrix A = Jt * J; // 3 x 3
        Matrix b = Jt * R; // 3 x 1

        Matrix delta = solve(A, b); // 3 x 1
        Vec3d delta_vec = mat2vec(delta);

        params = params - delta_vec * 0.01;
    }    
}