#include "model.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

const std::vector<std::vector<double>> alpha =
{
    {1.0/5.0},
    {3.0/40.0, 9.0/40.0},
    {44.0/45.0, -56.0/15.0, 32.0/9.0},
    {19372.0/6561.0, -25360.0/2187.0,	64448.0/6561.0,	-212.0/729.0},
    {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0},
    {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, 2187.0/6784.0, 11.0/84.0}
};

const std::vector<double> b = 
{
    35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0
};

std::vector<Vec3d> compute_accelerations(const std::vector<Object>& objects)
{
    std::vector<Vec3d> accelerations(objects.size());

    for (int i = 0; i < objects.size(); ++i)
    {
        for (int j = 0; j < objects.size(); ++j)
        {
            if (i == j) continue;

            Vec3d radius    = (objects[j].position - objects[i].position);
            double distance_sq = radius * radius;
            double distance = std::sqrt(distance_sq);
            if (distance < 1e-10) continue;

            Vec3d direction = radius / distance;
            accelerations[i] += G * objects[j].mass / distance_sq * direction;
        }
    }

    return accelerations;
}

Matrix compute_change_rate(const std::vector<Object>& objects, const Matrix& change_rate)
{
    Matrix res = compute_gravitational_gradient(objects);

    res = res * change_rate;

    return res;
}

Matrix compute_gravitational_gradient(const std::vector<Object>& objects)
{
    Matrix res(3);
    res.zeros();

    Object oumuamua = objects[0];

    for (int i = 1; i < objects.size(); ++i)
    {
        Object cur_obj = objects[i];

        Vec3d radius = oumuamua.position - cur_obj.position;
        double distance_sq = radius * radius;
        double dist = sqrt(distance_sq);
        
        double general = (G * cur_obj.mass) / pow(dist, 5);

        res.mat[0][0] += general * (-2 * pow(radius.x, 2) + pow(radius.y, 2) + pow(radius.z, 2));
        res.mat[0][1] += general * (-3 * radius.x * radius.y);
        res.mat[0][2] += general * (-3 * radius.x * radius.z);

        // res.mat[1][0] += general * (-3 * radius.x * radius.y);
        res.mat[1][1] += general * (pow(radius.x, 2) -2 * pow(radius.y, 2) + pow(radius.z, 2));
        res.mat[1][2] += general * (-3 * radius.z * radius.y);

        // res.mat[2][0] += general * (-3 * radius.x * radius.z);
        // res.mat[2][1] += general * (-3 * radius.z * radius.y);
        res.mat[2][2] += general * (pow(radius.x, 2) + pow(radius.y, 2) -2 * pow(radius.z, 2));
    }
    
    res.mat[1][0] = res.mat[0][1];
    res.mat[2][0] = res.mat[0][2];
    res.mat[2][1] = res.mat[1][2];

    return res;
}

SystemState derivative(const SystemState& state, const std::vector<Object>& objects)
{
    SystemState deriv;

    deriv.positions = state.velocities;

    deriv.velocities  = compute_accelerations(objects);

    deriv.change_rate = compute_change_rate(objects, state.change_rate);

    return deriv;
}

std::vector<Object> dopri5(std::vector<Object> objects, SystemState& state, double dt=1e-3)
{
    SystemState initial_state = state;

    SystemState k1 = derivative(initial_state, objects);


    SystemState temp = initial_state + k1 * alpha[0][0] * dt;
    SystemState k2 = derivative(temp, objects);

    temp = initial_state + (
        k1 * alpha[1][0] + 
        k2 * alpha[1][1]) * dt;
    SystemState k3 = derivative(temp, objects);

    temp = initial_state + (
        k1 * alpha[2][0] + 
        k2 * alpha[2][1] + 
        k3 * alpha[2][2]) * dt;
    SystemState k4 = derivative(temp, objects);

    temp = initial_state + (
        k1 * alpha[3][0] + 
        k2 * alpha[3][1] + 
        k3 * alpha[3][2] + 
        k4 * alpha[3][3]) * dt;
    SystemState k5 = derivative(temp, objects);

    temp = initial_state + (
        k1 * alpha[4][0] + 
        k2 * alpha[4][1] + 
        k3 * alpha[4][2] + 
        k4 * alpha[4][3] + 
        k5 * alpha[4][4]) * dt;
    SystemState k6 = derivative(temp, objects);

    temp = initial_state + (
        k1 * alpha[5][0] + 
        k2 * alpha[5][1] + 
        k3 * alpha[5][2] + 
        k4 * alpha[5][3] + 
        k5 * alpha[5][4] + 
        k6 * alpha[5][5]) * dt;
    SystemState k7 = derivative(temp, objects);

    SystemState result;
    result.positions.resize(objects.size());
    result.velocities.resize(objects.size());

    for (int i = 0; i < objects.size(); ++i)
    {
        result.positions[i] = initial_state.positions[i] + (
            b[0] * k1.positions[i] +
            b[1] * k2.positions[i] +
            b[2] * k3.positions[i] +
            b[3] * k4.positions[i] +
            b[4] * k5.positions[i] +
            b[5] * k6.positions[i] +
            b[6] * k7.positions[i]
        ) * dt;

        result.velocities[i] = initial_state.velocities[i] + (
            b[0] * k1.velocities[i] +
            b[1] * k2.velocities[i] +
            b[2] * k3.velocities[i] +
            b[3] * k4.velocities[i] +
            b[4] * k5.velocities[i] +
            b[5] * k6.velocities[i] +
            b[6] * k7.velocities[i]
        ) * dt;
    }

    result.change_rate = initial_state.change_rate + (
                        k1.change_rate * b[0] +
                        k2.change_rate * b[1] +
                        k3.change_rate * b[2] +
                        k4.change_rate * b[3] +
                        k5.change_rate * b[4] +
                        k6.change_rate * b[5] +
                        k7.change_rate * b[6] 
                    ) * dt;

    state = result;
    
    for (int i = 0; i < objects.size(); ++i)
    {
        objects[i].position = result.positions[i];
        objects[i].velocity = result.velocities[i];
    }

    return objects;
}

void init_state(std::vector<Object> objects, SystemState& state, const Matrix& change_rate_init)
{
    state.change_rate = change_rate_init;

    for (const auto& obj: objects)
    {
        state.positions.push_back(obj.position);
        state.velocities.push_back(obj.velocity);
    }
}

void write_state(const SystemState& state, std::ofstream& file)
{
    file << std::setprecision(15);
    for (const auto& pos: state.positions)
    {
        file << pos.x << " " << pos.y << " " << pos.z << "|";
    }

    file << "\n";
}

void integrate(
    std::vector<Object> objects,
    std::vector<SystemState>& states,
    std::vector<std::vector<Object>>& objects_trajectories,
    const Matrix& change_rate_init, 
    double t, double dt)
{   
    std::ofstream file("trajectory.txt");
    SystemState current_state;
    init_state(objects, current_state, change_rate_init);

    states.push_back(current_state);
    objects_trajectories.push_back(objects);

    for (double time = 0.0; time <= t; time += dt)
    {
        write_state(current_state, file);
        objects = dopri5(objects, current_state, dt);
        current_state.time = time + dt;

        objects_trajectories.push_back(objects);
        states.push_back(current_state);
    }
}
