#ifndef SYSTEMOFPARTICLES_H
#define SYSTEMOFPARTICLES_H
class Force;
class Timer;

class SystemOfParticles
{

    private:
    	int number_of_particles;
        double *r, *v, *F, *m;  // position, velocity, force, mass
        double T;               // Temperature
        double density;
        double volume;
        double dt;				// Time step
        Force *f;
        std::string *particle_name;
        
        
        double lateral_size;
        double lattice_parameter;
        
        void compute_interations();
        void move_particles();
        void compute_velocities();
        void check_wall_collisions();
        
        double knetic_energy();
        double potential_energy();
        
        void set_positions();
        void set_velocities(double mean = 0.0, double dispertion = 1.0);
        void set_equal_masses(double mass = 1.0);
        
        void set_particles_names();
        
        void show_infos(int, int, int);
        void store_files(int, int, int,
						std::string k_energy_file =  "kEnergy.dat",
						std::string p_energy_file =  "pEnergy.dat",
						std::string e_energy_file =  "Energy.dat");
		void store_xyz_file(bool append = false,
							std::string file_name = "positions.xyz",
							std::string delimiter = "\t");    
    public:
    	SystemOfParticles ();
        SystemOfParticles (int, double temperature = 300.0, double dens = 40.0, double time_step = 1.0e-4); 
        
        void set_initial_state(double mass = 1.0, double mean = 0.0, double dispertion = 1.0);
        
        void execute_interations(int, Timer *timer = nullptr);
        
        void load_masses(std::string);
        
        void load_names(std::string file = "particle_names.dat");
        
        void load_state(std::string);
        
        void store_state(std::string);
        
};

#endif
