#include <iostream>
#include "model.hpp"
#include "utils.hpp"


std::vector<Object> create_model_system() 
{
    std::vector<Object> objects;
    
    const double SUN_MASS = 1.989e30;
    const double MERCURY_MASS = 3.301e23;
    const double VENUS_MASS = 4.867e24;
    const double EARTH_MASS = 5.972e24;
    const double MARS_MASS = 6.417e23;
    const double JUPITER_MASS = 1.898e27;
    const double SATURN_MASS = 5.683e26;
    const double URANUS_MASS = 8.681e25;
    const double NEPTUNE_MASS = 1.024e26;

    Object oumuamua, sun, jupiter, earth;

    oumuamua.position = Vec3d(
        4.576215229884559E+09, 
        7.556775107045462E+08,    
        1.897489304528461E+09
    ) * 1000;
    oumuamua.velocity = Vec3d(
        2.486030535763478E+01, 
        3.705938998982580E+00,
        1.086798534442793E+01
    ) * 1000;
    oumuamua.mass = 4 * 10e11;

    sun.position = Vec3d(
        -1.351419336613188E+06,
        -1.263853537870709E+04,
        3.158480284404224E+04
    ) * 1000;
    sun.velocity = Vec3d(
        2.076458379360649E-03,
        -1.548435725926448E-02,
        7.866164391471479E-05
    ) * 1000;
    sun.mass = SUN_MASS;

    jupiter.position = Vec3d(
        -1.351419336613188E+06,
        -1.263853537870709E+04,
        3.158480284404224E+04
    ) * 1000;
    jupiter.velocity = Vec3d(
        2.076458379360649E-03,
        -1.548435725926448E-02,
        7.866164391471479E-05
    ) * 1000;
    jupiter.mass = JUPITER_MASS;

    earth.position = Vec3d(
        -1.351419336613188E+06,
        -1.263853537870709E+04, 
        3.158480284404224E+04
    ) * 1000;
    earth.velocity = Vec3d(
        2.076458379360649E-03,
        -1.548435725926448E-02,
        7.866164391471479E-05
    ) * 1000;
    earth.mass = EARTH_MASS;

    objects.push_back(oumuamua);
    objects.push_back(sun);
    objects.push_back(jupiter);
    objects.push_back(earth);

    return objects;
}


int main(int argc, char* argv[])
{
    std::vector<Object> system = create_model_system();


    return 0;
}