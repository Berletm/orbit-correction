#include "correction.hpp"
#include "model.hpp"
#include "interpolation.hpp"
#include "gauss.hpp"
#include <cmath>
#include <algorithm>

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

Matrix calculate_state_change(const Vec3d& computed)
{
    Matrix res(2, 3);

    double general = pow(computed.x, 2) + pow(computed.y, 2);
    double dist_sq = computed * computed;
    res.mat[0][0] = - computed.y / general;
    res.mat[0][1] = computed.x / general;
    res.mat[0][2] = 0;

    res.mat[1][0] = - computed.x * computed.z / (sqrt(general) * dist_sq);
    res.mat[1][1] = - computed.y * computed.z / (sqrt(general) * dist_sq);
    res.mat[1][2] = sqrt(general) / dist_sq;

    return res;
}

Matrix calculate_jacobian(const Matrix& state_change, const Matrix& change_rate)
{
    Matrix j = state_change * change_rate;

    return j;
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

        res.mat[0 + shift][0] = cur_vec.ra; 
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

void correction(std::vector<Object>& initial_state, const std::vector<Celestial>& observed, const std::vector<Vec3d>& obs_position, const std::vector<double>& obs_time)
{   
    Vec3d& params = initial_state[0].position;

    double t_end = *std::max_element(obs_time.begin(), obs_time.end());

    std::vector<SystemState> states;
    std::vector<std::vector<Object>> objects_trajectories;
    Matrix change_rate_init(3);
    change_rate_init.identity();
    integrate(initial_state, states, objects_trajectories, change_rate_init, t_end, 1e-3);

    std::vector<SystemState> states_at_measurements;

    for (const auto& time: obs_time)
    {
        SystemState state = interpolate(objects_trajectories, states, time);
        states_at_measurements.push_back(state);
    }

    std::vector<Vec3d> computed;
    for (int i = 0; i < states_at_measurements.size(); ++i)
    {
        SystemState state = states_at_measurements[i];
        Vec3d observatory_vec = obs_position[i]; 

        Vec3d oumuamua_pos = state.positions[0];
        Vec3d earth_pos    = state.positions[1];

        Vec3d temp_vec     =  earth_pos + observatory_vec;

        Vec3d computed_vec = oumuamua_pos - temp_vec; // not sure about operands order 
        computed.push_back(computed_vec);
    }

    std::vector<Celestial> computed_angles = cart2celestial(computed);

    std::vector<Celestial> r = calculate_residuals(observed, computed_angles);

    std::vector<Matrix> jacobians;

    for (int i = 0; i < states_at_measurements.size(); ++i)
    {
        SystemState state = states_at_measurements[i];
        Vec3d comp = computed[i];
        Matrix state_change = calculate_state_change(comp);
        Matrix j = calculate_jacobian(state_change, state.change_rate);

        jacobians.push_back(j);
    }

    Matrix J = stack_matrix(jacobians); // 2N x 3 (rows x cols)
    Matrix R = stack_vector(r);         // 2N x 1 (rows x cols)

    Matrix Jt = J; // 3 x 2N
    Jt.transpose();

    Matrix A = Jt * J; // 3 x 3
    Matrix b = Jt * R; // 3 x 1

    b = b * -1;

    Matrix delta = solve(A, b); // 3 x 1

    Vec3d delta_vec = mat2vec(delta);

    params = params - delta_vec;
}