#include <iostream>
#include "model.hpp"
#include "utils.hpp"
#include "correction.hpp"


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
        1.452674920249601E+08, 
        7.476202044512266E+07,    
        -1.071281870832360E+07
    ) * 1000;
    oumuamua.velocity = Vec3d(
        4.483537051061160E+01,
        1.039892688342698E+01,
        1.433813590587155E+01
    ) * 1000;
    oumuamua.mass = 1e11;

    sun.position = Vec3d(
        3.348986030140055E+05,
        8.546522676943250E+05,
        -1.952911753764993E+04
    ) * 1000;
    sun.velocity = Vec3d(
        -8.992512340765906E-03,
        9.545307638330358E-03,
        2.103021639645862E-04
    ) * 1000;
    sun.mass = SUN_MASS;

    jupiter.position = Vec3d(
        -6.876378394323186E+08,
        -4.342916936111128E+08,
        1.718173022009712E+07
    ) * 1000;
    jupiter.velocity = Vec3d(
        6.823103108737919E+00,
        -1.042433580371820E+01,
        -1.092815489926169E-01
    ) * 1000;
    jupiter.mass = JUPITER_MASS;

    earth.position = Vec3d(
        1.400225151174905E+08,
        5.334382691880187E+07,
        -2.224282551313937E+04
    ) * 1000;
    earth.velocity = Vec3d(
        -1.096275218788954E+01,
        2.779127096624116E+01,
        -1.998546477564034E-03
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

    std::vector<Vec3d> observers;
    std::vector<Celestial> observed;
    std::vector<double> time;

    read_observed_data(time, observers, observed);
    
    correction(system, observed, observers, time);

    // std::vector<SystemState> states;
    // std::vector<std::vector<Object>> objects_trajectories;
    // Matrix change_rate_init(3);
    // change_rate_init.identity();
    // system[0].position = Vec3d(-1.29024e+55, -3.66654e+54, 2.28563e+53);
    // integrate(system, states, objects_trajectories, change_rate_init, 365*24*60*60, 60*60*24); // step = 6 hours

    return 0;
}