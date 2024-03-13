cmake_minimum_required(VERSION 3.14)
project(moldyn CXX)

option(ENABLE_OPENMP "use OpenMP for parallel" OFF)

set(HEADERS headers/fileHandler.hpp headers/Force.hpp headers/SystemOfParticles.hpp headers/Timer.hpp headers/eigen_wrapper.hpp)
set(SOURCES sources/Force.cpp sources/SystemOfParticles.cpp sources/Timer.cpp main.cpp)

add_executable(moldyn ${SOURCES} ${HEADERS})

target_include_directories(moldyn PUBLIC headers)
install(TARGETS moldyn)

if (ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(moldyn PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(STATUS "OpenMP not found!")
    endif()
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(moldyn PRIVATE Eigen3::Eigen)
