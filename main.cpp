#include <iostream>
#include "SystemOfParticles.hpp"
#include "etime.h"


int main() {
    
    SystemOfParticles gas(27, 180.0);
    gas.set_initial_state(5.0);

    etime timer;
    
    timer.start();
    
    	gas.execute_interations(1000000);
    	
    timer.end("Simulation lasted ");
    
    return 0;
}
