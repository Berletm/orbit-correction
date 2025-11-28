#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <vector>

struct Vec3d
{
    double x, y, z;

    inline Vec3d(): x(0), y(0), z(0) {}

    inline Vec3d(double x, double y, double z): x(x), y(y), z(z) {}

    inline Vec3d operator+(const Vec3d& other) const
    {
        return Vec3d(
            this->x + other.x, 
            this->y + other.y, 
            this->z + other.z);
    }

    inline void operator+=(const Vec3d& other)
    {
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
    }

    inline Vec3d operator-(const Vec3d& other) const
    {
        return Vec3d(
            this->x - other.x, 
            this->y - other.y, 
            this->z - other.z);
    }

    inline double operator*(const Vec3d& other) const
    {
        return this->x * other.x + this->y * other.y + this->z * other.z;
    }

    inline Vec3d operator*(const double scalar) const
    {
        return Vec3d(this->x * scalar, this->y * scalar, this->z * scalar);
    }

    inline Vec3d operator/(const double scalar) const
    {
        return Vec3d(
            this->x / scalar, 
            this->y / scalar, 
            this->z / scalar);
    }

    friend Vec3d operator*(double scalar, const Vec3d& vec);
};

inline Vec3d operator*(double scalar, const Vec3d& vec) 
{
    return Vec3d(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

struct SystemState
{
    std::vector<Vec3d> positions;
    std::vector<Vec3d> velocities;

    SystemState operator+(const SystemState& other)
    {
        SystemState res;
        res.positions.resize(positions.size());
        res.velocities.resize(velocities.size());

        for (int i = 0; i < positions.size(); ++i)
        {
            res.positions[i] = this->positions[i] + other.positions[i];
            res.velocities[i] = this->velocities[i] + other.velocities[i];
        }

        return res;
    }

    SystemState operator*(double scalar)
    {
        SystemState res;
        res.positions.resize(positions.size());
        res.velocities.resize(velocities.size());

        for (int i = 0; i < positions.size(); ++i)
        {
            res.positions[i] = res.positions[i] * scalar;
            res.velocities[i] = res.velocities[i] * scalar;
        }

        return res;
    }

};

struct Object
{
    Vec3d position;
    Vec3d velocity;
    Vec3d acceleration;
    
    double mass;
};


#endif // UTILS_HPP
