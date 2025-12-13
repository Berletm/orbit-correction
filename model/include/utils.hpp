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

struct Celestial
{
    double ra, dec;
    inline Celestial(double ra, double dec): ra(ra), dec(dec) {}

    inline Celestial operator-(const Celestial& other) const
    {
        return Celestial(this->ra - other.ra, this->dec - other.dec);
    }
};

struct Matrix
{
    std::vector<std::vector<double>> mat;
    int cols, rows;

    Matrix(int cols, int rows): cols(cols), rows(rows)
    {
        mat.resize(rows);

        for (auto& row: mat)
        {
            row.resize(cols, 0);
        }
    }

    Matrix(int size): cols(size), rows(size)
    {
        Matrix(size, size);
    }

    inline void identity()
    {
        for (int i = 0; i < rows; ++i)
        {
            mat[i][i] = 1;
        }
    }

    inline void zeros()
    {
        for (auto& row: mat)
        {
            std::fill(row.begin(), row.end(), 0);
        }
    }

    inline Matrix operator*(const Matrix& other) const
    {   // (n * k) * (k * p) = (n * p)
        Matrix res(this->rows, other.cols);

        for (int i = 0; i < this->rows; ++i)
        {
            for (int j = 0; j < other.cols; ++j)
            {
                for (int k = 0; k < other.rows; ++k)
                {
                    res.mat[i][j] += this->mat[i][k] * other.mat[k][j];
                }
            }
        }
    }

};

struct SystemState
{
    std::vector<Vec3d> positions;
    std::vector<Vec3d> velocities;
    Matrix change_rate;

    SystemState(): change_rate(3) { change_rate.identity(); }

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
