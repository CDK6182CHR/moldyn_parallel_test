#include <iostream>
#include "etime.hpp"
#include "SystemOfParticles.hpp"



int main() {
    
    std::cout << "\n\t ********************\n";
    std::cout << "\n\t Program Moldyn V 1.1\n";
    std::cout << "\n\t ********************\n";
    
    std::cout << "\n\t A simple Particle Dynamics program \n";
    std::cout << "\t for simulating interacting point particles.\n";
    
    std::cout << "\n\tCopyright (C) 2019  Vagner Bessa, Crateus - CE, Brazil.\n";
    
    SystemOfParticles gas(64, 80.0);
    gas.set_initial_state(3.9);

    etime timer;
    
    timer.start();
    
    	gas.execute_interations(1000000);
    	
    timer.end("Simulation lasted ");
    
    return 0;
}
