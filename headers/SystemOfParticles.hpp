#ifndef SYSTEMOFPARTICLES_H
#define SYSTEMOFPARTICLES_H

#include <memory>
#include "Force.hpp"

class Force;
class Timer;

class SystemOfParticles
{

private:
	int number_of_particles;
	std::unique_ptr<double[]>  r, v, F, m;  // position, velocity, force, mass
	double T;               // Temperature
	double density;
	double volume;
	double dt;				// Time step
	double pressure;
	double K_energy, P_energy, energy;
	std::unique_ptr<Force> f;
	std::unique_ptr<std::string[]> particle_name;


	double lateral_size;
	double lattice_parameter;

	void compute_interations();
	void move_particles();
	void compute_velocities();
	double check_wall_collisions();   // return the pressure due to collisions

	double knetic_energy();
	double Temperature();
	double potential_energy();

	void set_positions();
	void set_velocities(double mean = 0.0, double dispertion = 1.0);
	void set_equal_masses(double mass = 1.0);

	void show_infos(int, int, int, double, double);
	void store_files(int, int, int, double, double,
		std::string k_energy_file = "kEnergy.dat",
		std::string p_energy_file = "pEnergy.dat",
		std::string e_energy_file = "Energy.dat",
		std::string temperature_file = "Temperature.dat",
		std::string average_pressure_file = "Avg_Pressure.dat",
		std::string heat_capacity_file = "heat_capacity.dat");
	void store_xyz_file(bool append = false,
		std::string file_name = "positions.xyz",
		std::string delimiter = "\t");
public:

	SystemOfParticles();

	SystemOfParticles(int, double temperature = 300.0, double dens = 40.0, double time_step = 1.0e-4);

	void set_initial_state(double mass = 1.0, double mean = 0.0, double dispertion = 1.0);

	void set_particles_name(std::string name = "Ar");

	void execute_interations(int, Timer* timer = nullptr);

	void load_state(std::string);

	void store_state(std::string);

};

#endif
