#include <iostream>
#include "model.hpp"
#include "utils.hpp"


std::vector<Object> create_solar_system_precise() 
{
    std::vector<Object> objects;
    
    double sun_mass = 1.989e30;
    double earth_mass = 5.972e24;
    
    Object sun;
    sun.mass = sun_mass;
    sun.position = Vec3d(0, 0, 0);
    sun.velocity = Vec3d(0, 0, 0);
    sun.acceleration = Vec3d(0, 0, 0);
    
    Object earth;
    earth.mass = earth_mass;
    
    earth.position = Vec3d(
        1.466736208258275e11, 
        -2.573616899222266e10,
        1.478616600000000e7   
    );
    
    earth.velocity = Vec3d(
        4.982246315320939e3,  
        2.927048310824395e4,  
        -3.346551755298911e-1  
    );
    
    earth.acceleration = Vec3d(0, 0, 0);
    
    objects.push_back(sun);
    objects.push_back(earth);
    
    return objects;
}


int main(int argc, char* argv[])
{
    std::vector<Object> system = create_solar_system_precise();


    double total_time = 365.25 * 24 * 3600;
    double dt = 3600;

    integrate(system, total_time, dt);

    return 0;
}