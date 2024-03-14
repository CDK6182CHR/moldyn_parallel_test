cmake_minimum_required(VERSION 3.22.3)
project(moldyn CXX)

option(ENABLE_OPENMP "use OpenMP for parallel" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
    add_compile_options(/Zc:__cplusplus)
endif ()


set(HEADERS headers/fileHandler.hpp headers/Force.hpp headers/SystemOfParticles.hpp headers/Timer.hpp)
set(SOURCES sources/Force.cpp sources/SystemOfParticles.cpp sources/Timer.cpp main.cpp)



add_executable(moldyn ${SOURCES} ${HEADERS})
install(TARGETS moldyn)

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

find_package(IntelDPCPP REQUIRED)
message(STATUS "Found SYCL include path: ${SYCL_INCLUDE_DIR}")


