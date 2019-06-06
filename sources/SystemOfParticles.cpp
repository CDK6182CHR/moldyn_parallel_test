#include <time.h>
#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "etime.hpp"
#include "Force.hpp"
#include "SystemOfParticles.hpp"
#include "fileHandler.hpp"



SystemOfParticles::SystemOfParticles () : 
	r(nullptr),
	v(nullptr),
	m(nullptr),
	number_of_particles (0),
	f(nullptr),
	particle_name(nullptr),
	dt(0.0),
	lateral_size(0.0),
	lattice_parameter(0.0)
	{}
        
SystemOfParticles::SystemOfParticles (int N, double temperature, double time_step) {
	
	number_of_particles = N;
	r = new double[3*N];	//Position
	v = new double[3*N];	//Velocity
	F = new double[3*N];	//Force
	m = new double[N];	//Mass
	particle_name = new std::string [N];
	f = new Force();		//Force of interaction
	T = temperature/f->unit_of_temperature();
	dt = time_step;
}

void SystemOfParticles::set_initial_state(double L, double mass, double mean, double dispertion) {
	lateral_size = L;

	set_particles_names();

	set_equal_masses(mass);
	set_positions();
	set_velocities(mean, dispertion);

	compute_interations(); // In order to aquire units;
	for (unsigned int i = 0; i < 3*number_of_particles; i += 1) {
	    F[i] = 0.0;
	}
	
	
	std::cout << "\n\tSystem of units used:\n";
	std::cout << "\tUnit of energy: " << f->unit_of_energy() << std::endl;
	std::cout << "\tUnit of space: " << f->unit_of_space() << std::endl;
	std::cout << "\tUnit of temperature: " << f->unit_of_temperature() << std::endl;
	std::cout << "\tUnit of time: " << f->unit_of_time() << std::endl;
	
}

void SystemOfParticles::execute_interations(int number_of_interations, etime *timer) {
	
	int factor_percent = number_of_interations/100;
	int factor_ecran = factor_percent;
	int factor_store_state = 100;
	int factor_xy = 1;
	
	if (timer == nullptr) {
	    timer = new etime();
	}
	
	timer->start();
	for (unsigned int it = 0; it < number_of_interations; it += 1) {
		
		show_infos(it, factor_ecran, factor_percent);
		
		for (unsigned int i = 0; i < 3*number_of_particles; i += 1) {
		    F[i] = 0.0;
		}	
		
		//================================================
		//      THE VELOCITY-VERLET BLOCK
		compute_interations();

		move_particles();
		
		compute_interations();
		      
		compute_velocities();
		//================================================
		
		check_wall_collisions();
		
		store_files(it, factor_store_state, factor_xy);
		
		if (!(it%factor_ecran)) timer->register_time("simulation time: ");
		
	}
	timer->end("simulation time: ");

}

void SystemOfParticles::compute_interations() {
	
	// computes N(N - 1) pairs of interactions 
	// the force array is not set to zero in order to compute
	// the velocity at the last step of the Velocity-Verlet
	// algorithm
	double *F_ptr;
    for (unsigned int i = 0; i < 3*number_of_particles - 3; i += 3) {
        for (unsigned int j = i + 3; j < 3*number_of_particles; j += 3) {
        
            F_ptr = (*f)(&r[i], &r[j], m[i/3], m[j/3]);
            for (unsigned int k = 0; k < 3; k += 1) {
            
                F[i + k] += F_ptr[k];
                F[j + k] -= F_ptr[k];
                
            }
            
        }
    }
}

void SystemOfParticles::move_particles() {

    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
    	for (unsigned int j = 0; j < 3; j += 1) {
    	
    		r[i + j] += v[i + j]*dt + 0.5*(F[i + j]/m[i/3])*dt*dt;
    		
    	}
   	}
}

void SystemOfParticles::compute_velocities() {

    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
    	for (unsigned int j = 0; j < 3; j += 1) {
    	
    		v[i + j] += 0.5*(F[i + j]/m[i/3])*dt;
    		
    	}
	}
}

double SystemOfParticles::knetic_energy() {
    double v2 = 0.0;
    
    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
        for (unsigned int j = 0; j < 3; j += 1) {
            v2 += m[i/3]*v[i + j]*v[i + j];
        }
    }
    
    return 0.5*v2;
    
}

double SystemOfParticles::potential_energy() {
    double U = 0.0;
    for (unsigned int i = 0; i < 3*number_of_particles - 3; i += 3) {
        for (unsigned int j = i + 3; j < 3*number_of_particles; j += 3) {
            U += f->potential(&r[i], &r[j], m[i/3], m[j/3]);
        }
    }
    return U;
}

void SystemOfParticles::set_positions() {

    int n;
    
    n = int(ceil(pow(number_of_particles,1.0/3.0))); // Number of molecules per direction
    lattice_parameter = lateral_size/n;        // Distance between molecules
    
    int idx = 0;
    int TN = 3*number_of_particles;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++){
            for (unsigned int k = 0; k < n; k++) {
                if (idx < TN){
                    r[idx] = (i + 0.5)*lattice_parameter;
                    r[idx + 1] = (j + 0.5)*lattice_parameter;
                    r[idx + 2] = (k + 0.5)*lattice_parameter;
                }
                idx += 3;
            }
        }
    }

}

void SystemOfParticles::set_velocities(double mean, double dispertion) {

	double Vcm[3] = {0.,0.,0.};
	double M = 0.0;
	double v2 = 0.0;
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::normal_distribution<double> distribution (mean, dispertion);
    
    for (unsigned int i = 0; i < 3*number_of_particles; i += 1) {
        v[i] = distribution(generator);
    }
    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
        for (unsigned int j = 0; j < 3; j += 1) {
            Vcm[j] += m[i/3]*v[i + j];
        }
    }
    
    for (unsigned int i = 0; i < number_of_particles; i += 1) {
        M += m[i];
    }
    
    for (unsigned int j = 0; j < 3; j += 1) {
        Vcm[j] /= number_of_particles*M;
    }

    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
        for (unsigned int j = 0; j < 3; j += 1) {
            v[i + j] -= Vcm[j];
            v2 += m[i/3]*v[i + j]*v[i + j];
        }
    }
    
    double To = 0.5*v2/(3*(number_of_particles - 1));
    double temperature_factor = sqrt(T/To);
    
    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
        for (unsigned int j = 0; j < 3; j += 1) {
            v[i + j] *= temperature_factor;
        }
    }   
    
}

void SystemOfParticles::set_equal_masses(double mass) {
    for (unsigned int i = 0; i < number_of_particles; i += 1) {
        m[i] = mass;
    }
}

void SystemOfParticles::set_particles_names() {
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		particle_name[i] = "He"; // Default
	}
}

void SystemOfParticles::check_wall_collisions() {
  // Elastic walls
    for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
        for (unsigned int j = 0; j < 3; j += 1) {
			if (r[i + j] < 0.0) {
				v[i + j] *=-1.; //- elastic walls
			}
			if (r[i + j] >= lateral_size) {
				v[i + j]*=-1.;  //- elastic walls
			}
		}
	}
}

void SystemOfParticles::load_masses(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		file >> m[i];
	}
	
	file.close();
}

void SystemOfParticles::load_names(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		file >> particle_name[i];
	}
	
	file.close();
}

void SystemOfParticles::load_state(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	file >> number_of_particles;
	file >> dt;
	
	int N = number_of_particles;
	
	r = new double[3*N];	//Position
	v = new double[3*N];	//Velocity
	F = new double[3*N];	//Force
	m = new double[N];	//Mass
	particle_name = new std::string [N];
	f = new Force();		//Force of interaction
		
	for (unsigned int i = 0; i < 3*number_of_particles; i += 3) {
		file >> particle_name[i/3] >> m[i/3] >> r[i] >> r[i + 1] >> r[i + 2] >> v[i] >> v[i + 1] >> v[i + 2];
	}
	
	file.close();
}

void SystemOfParticles::store_state(std::string file_name) {
	std::ofstream file;
	file.open(file_name);
	file << number_of_particles << std::endl;
	file << dt << std::endl;
	
	int N = number_of_particles;
	
	r = new double[3*N];	//Position
	v = new double[3*N];	//Velocity
	F = new double[3*N];	//Force
	m = new double[N];	//Mass
	particle_name = new std::string [N];
	f = new Force();		//Force of interaction
		
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		file << particle_name[i] << "\t";
		file << m[i]  << "\t";
		file << r[i] << "\t" << r[i + 1] << "\t" << r[i + 2] << "\t";
		file << v[i] << "\t" << v[i + 1] << "\t" << v[i + 2] << std::endl;
	}
	file.close();
}

void SystemOfParticles::store_xyz_file(bool append, std::string file_name, std::string delimiter) {
    std::ofstream file;
    if (append) {
		file.open(file_name, std::ios::app);
    } else {
    	file.open(file_name);
    	file << number_of_particles << std::endl;
    }
    for (unsigned int idx = 0; idx < 3*number_of_particles; idx += 3) {
    	int i = idx/3;
        file << particle_name[i];
        file << delimiter << r[idx] << delimiter << r[idx + 1] << delimiter << r[idx + 2] << std::endl; 
    }
    file.close();
}


void SystemOfParticles::store_files(int iterator, int factor_xyz, int factor_xy,
									std::string k_energy_file,
									std::string p_energy_file,
									std::string e_energy_file) {
	
	if (iterator == 0) {
		storeData(f->unit_of_energy(), k_energy_file);
		storeData(f->unit_of_energy(), p_energy_file);
		storeData(f->unit_of_energy(), e_energy_file);
		store_xyz_file();
	}
	
	if (!(iterator%factor_xyz)) {
		store_xyz_file(true);
	}
	if (!(iterator%factor_xy)) {
		double K = knetic_energy()*f->unit_of_energy();
		double P = potential_energy()*f->unit_of_energy();
		double E = K + P;
		double time = dt*iterator;
		storeXYData(time, K, k_energy_file);
		storeXYData(time, P, p_energy_file);
		storeXYData(time, E, e_energy_file);
	}
}
void SystemOfParticles::show_infos(int iterator, int factor_ecran, int factor_percent) {
	
	double K, P, E;
	
	std::cout << std::setprecision(1);
	
    if (!iterator) {
        std::cout << "\n\t========================================\n";
        std::cout << "\t|Simulation state \t| knetic energy \t| potential energy \t| Total energy \t| Temperature (K) \t\n\n \t\n\n";
        std::cout.flush();
    }
    if (!(iterator%factor_ecran)) {
//            if (!((it/fac_percent)%10)) {
            K = knetic_energy();
            P = potential_energy();
            E = K + P;
            std::cout << std::fixed << "\t|" << double(iterator/factor_percent) << "% \t        ";
            std::cout << std::scientific;
            std::cout << "\t|" << K*f->unit_of_energy() << "        " ;
           	std::cout << "\t|" << P*f->unit_of_energy() << "        " ;
            std::cout << "\t|" << E*f->unit_of_energy() << "     ";
            std::cout << "\t|" << K/(3.0*(number_of_particles - 1.0))*f->unit_of_temperature() << "        ";
            std::cout << std::endl;
            std::cout.flush();
//            }
    }

}
