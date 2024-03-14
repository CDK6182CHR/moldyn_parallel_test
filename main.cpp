#include <iostream>
#include "Timer.hpp"
#include "SystemOfParticles.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

void set_inputs(int*, double*, double*, std::string *input_file_name = nullptr);

int main() {
#ifdef _OPENMP
    std::cout << "Info: OpenMP is enabled with " << omp_get_max_threads() << " threads" << std::endl;
#endif
    std::cout << "__cplusplus: " << __cplusplus << std::endl;
    int number_of_particles = 216;
    double temperature = 80.0;
    double density = 3.5e4;
    
    set_inputs(&number_of_particles, &temperature, &density);
    
    SystemOfParticles gas(number_of_particles, temperature, density);
    gas.set_initial_state();

    Timer timer;
    
    timer.start_timer();
    
    	gas.execute_interations(20000);
    	
    timer.end_timer("Simulation lasted ");
    
    return 0;
}

void set_inputs(int* number_of_particles, double* temperature, double* density, std::string *input_file_name) {
    if (input_file_name == nullptr) {
        std::cout << "\n\t********************\n";
        std::cout << "\n\tProgram Moldyn V 1.1\n";
        std::cout << "\n\t********************\n";
        
        std::cout << "\n\tA simple Particle Dynamics program \n";
        std::cout << "\tfor simulating interacting point particles.\n";
        
        std::cout << "\n\tCopyright (C) 2019  Vagner Bessa, Crateus - CE, Brazil.\n";
        
        std::cout << "\n\tRuning with default values em parameters:\n";
        std::cout << "\t" << *number_of_particles << " atoms of Ar";
        std::cout << " at T = " <<  *temperature;
        std::cout << " and density of " << *density << " moles/m^3.\n";
    }
}
