#include "SystemOfParticles.hpp"

SystemOfParticles::SystemOfParticles () : r(nullptr), v(nullptr), m(nullptr) {}
        
SystemOfParticles::SystemOfParticles (int N) {
            r = new double[3*N];
            v = new double[3*N];
            m = new double[3*N];
}
