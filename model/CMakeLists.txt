cmake_minimum_required(VERSION 3.11...3.20)

project(orbit)

add_executable(orbit
    src/main.cpp    
    src/model.cpp
    src/correction.cpp
    src/interpolation.cpp
    src/gauss.cpp
)

target_include_directories(orbit PRIVATE include)
target_compile_options(orbit PRIVATE -fsanitize=address -g)
target_link_options(orbit PRIVATE -fsanitize=address -g)
