cmake_minimum_required(VERSION 3.11...3.20)

project(orbit)

add_executable(orbit
    src/main.cpp    
    src/model.cpp
)

target_include_directories(orbit PRIVATE include)
