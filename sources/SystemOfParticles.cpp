#include <time.h>
#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Timer.hpp"
#include "Force.hpp"
#include "SystemOfParticles.hpp"
#include "fileHandler.hpp"



SystemOfParticles::SystemOfParticles() :
	r(),
	v(),
	m(),
	number_of_particles(0),
	f(nullptr),
	particle_name(nullptr),
	dt(0.0),
	lateral_size(0.0),
	lattice_parameter(0.0)
{}

SystemOfParticles::SystemOfParticles(int N, double temperature, double dens, double time_step):
	r(N, 3), v(N, 3), F(N, 3), m(N)
{

	number_of_particles = N;
	particle_name.reset(new std::string[N]);
	f.reset(new Force());		            //Force of interaction

	T = temperature / f->unit_of_temperature();

	density = dens;
	volume = N / (density * f->get_avogadro_constant());
	std::cout << "\tVolume of the system: " << volume << " m³\n";
	volume /= pow(f->unit_of_space(), 3);
	lateral_size = pow(volume, 1.0 / 3.0);
	std::cout << "\tLateral size: " << lateral_size * f->unit_of_space() << " m\n";
	dt = time_step;

}

void SystemOfParticles::set_initial_state(double mass, double mean, double dispertion) {

	set_particles_name();

	set_equal_masses(mass);
	set_positions();
	set_velocities(mean, dispertion);

	compute_interations(); // In order to aquire units;
	F.fill(0);

	std::cout << "\n\tSystem of units used (in S.I. units):\n";
	std::cout << "\tUnit of energy: " << f->unit_of_energy() << std::endl;
	std::cout << "\tUnit of space: " << f->unit_of_space() << std::endl;
	std::cout << "\tUnit of time: " << f->unit_of_time() << std::endl;
	std::cout << "\tUnit of temperature: " << f->unit_of_temperature() << std::endl;


}

void SystemOfParticles::execute_interations(int number_of_interations, Timer* timer) {

	int factor_percent = number_of_interations / 100;
	int factor_ecran = factor_percent;
	int factor_store_state = 1000;
	int factor_xy = 1000;

	double average_temperature = 0.0;
	double average_pressure = 0.0;
	double average_energy = 0.0;
	double square_energy = 0.0;
	double heat_capacity = 0.0;
	double time;
	double gas_constant = 0.0;

	if (timer == nullptr) {
		timer = new Timer();
	}

	timer->start_timer();
	for (unsigned int it = 1; it < number_of_interations + 1; it += 1) {

		double P_loc = 0, K_loc = 0, U_loc = 0;
		F.fill(0);

#ifdef _OPENMP
#pragma omp parallel reduction(+:P_loc) reduction(+:K_loc) reduction(+:U_loc)
#endif
		{
//#ifdef _OPENMP
//#pragma omp for
//#endif
//			for (int i = 0; i < 3 * number_of_particles; i += 1) {
//				F[i] = 0.0;
//			}

			//================================================
			//      THE VELOCITY-VERLET BLOCK
			compute_interations();

			move_particles();

			compute_interations();

			compute_velocities();
			//================================================


			P_loc = check_wall_collisions();

			K_loc = knetic_energy();
			U_loc = potential_energy();
		}   // omp parallel

		pressure = P_loc;
		K_energy = K_loc;
		P_energy = U_loc;
		T = K_energy * 2 / (3 * (number_of_particles - 1));
		energy = K_energy + P_energy;

		average_pressure += pressure;
		average_temperature += T;
		average_energy += energy;
		square_energy += energy * energy;
		time = it * dt;
		heat_capacity = (square_energy / it - pow(average_energy / it, 2)) / (T * T);        // From canonical ensemble

		gas_constant = average_pressure * volume / (number_of_particles * average_temperature);
		gas_constant *= f->unit_of_gas_constant();

		store_files(it, factor_store_state, factor_xy, average_pressure, heat_capacity);

		show_infos(it, factor_ecran, factor_percent, average_pressure, gas_constant);

		if (!(it % factor_ecran)) timer->register_time("simulation time: ");

	}
	timer->end_timer();

}

void SystemOfParticles::compute_interations() {

	// computes N(N - 1) pairs of interactions 
	// the force array is not set to zero in order to compute
	// the velocity at the last step of the Velocity-Verlet
	// algorithm
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for (int i = 0; i < number_of_particles; i++) {
		for (int j = 0; j < number_of_particles; j++) {
			if (i == j) continue;

			auto F_sub = f->operator()(r.row(i), r.row(j));   // move construct
			F.row(i) += F_sub;   // assign

		}
	}
}

void SystemOfParticles::move_particles() {

	// Eigen vectorized operation
	r += v * dt + F.colwise() / m * dt * dt * 0.5;
}

void SystemOfParticles::compute_velocities() {

	// Eigen vectorized operation
	v += F.colwise() / m * 0.5 * dt;
}

double SystemOfParticles::check_wall_collisions() {
	// Elastic walls
	double P = 0.0;

#ifdef _OPENMP
#pragma omp for 
#endif
	for (int i = 0; i < number_of_particles; i++) {
		for (int j = 0; j < 3; j += 1) {
			if (r(i, j) < 0.0) {
				v(i, j) *= -1.; //- elastic walls
				P += 2.0 * m[i] * abs(v(i, j)) / dt;        // the average forçe from the walls
			}
			if (r(i, j) >= lateral_size) {
				v(i, j) *= -1.;  //- elastic walls
				P += 2.0 * m[i] * abs(v(i, j)) / dt;        // the average forçe from the walls
			}
		}
	}

	return P / (6.0 * lateral_size * lateral_size);            // the average forçe form the walls over area

}

double SystemOfParticles::knetic_energy() {
	double v2 = 0.0;

#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < number_of_particles; i++) {
		for (int j = 0; j < 3; j++) {
			v2 += m[i] * v(i, j) * v(i, j);
		}
	}
	//T = v2/(3*(number_of_particles - 1));   // 2024.03.09  for parallel version: put outside
	return 0.5 * v2;

}

double SystemOfParticles::Temperature() {
	T = 2.0 * knetic_energy() / (3 * number_of_particles);
	return T;
}

double SystemOfParticles::potential_energy() {
	double U = 0.0;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for (int i = 0; i < number_of_particles - 1; i++) {
		for (int j = i + 1; j < number_of_particles; j++) {
			U += f->potential(r.row(i), r.row(j));
		}
	}
	return U;
}

void SystemOfParticles::set_positions() {

	int n;

	n = int(ceil(pow(number_of_particles, 1.0 / 3.0))); // Number of molecules per direction
	lattice_parameter = lateral_size / n;        // Distance between molecules

	int idx = 0;
	int TN = number_of_particles;
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n; j++) {
			for (unsigned int k = 0; k < n; k++) {
				if (idx < TN) {
					r(idx, 0) = (i + 0.5) * lattice_parameter;
					r(idx, 1) = (j + 0.5) * lattice_parameter;
					r(idx, 2) = (k + 0.5) * lattice_parameter;
				}
				idx++;
			}
		}
	}

}

void SystemOfParticles::set_velocities(double mean, double dispertion) {

	double Vcm[3] = { 0.,0.,0. };
	double M = 0.0;
	double v2 = 0.0;
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(mean, dispertion);

	auto v_view = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
	for (unsigned int i = 0; i < 3 * number_of_particles; i++) {
		v_view[i] = distribution(generator);
	}
	for (unsigned int i = 0; i < number_of_particles; i++) {
		for (unsigned int j = 0; j < 3; j += 1) {
			Vcm[j] += m[i] * v(i, j);
		}
	}

	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		M += m[i];
	}

	for (unsigned int j = 0; j < 3; j += 1) {
		Vcm[j] /= number_of_particles * M;
	}

	for (unsigned int i = 0; i < number_of_particles; i++) {
		for (unsigned int j = 0; j < 3; j += 1) {
			v(i, j) -= Vcm[j];
			v2 += m[i] * v(i, j) * v(i, j);
		}
	}

	double To = v2 / (3 * (number_of_particles - 1));
	double temperature_factor = sqrt(T / To);

	for (unsigned int i = 0; i < number_of_particles; i ++) {
		for (unsigned int j = 0; j < 3; j += 1) {
			v(i, j) *= temperature_factor;
		}
	}

}

void SystemOfParticles::set_equal_masses(double mass) {
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		m[i] = mass;
	}
}

void SystemOfParticles::set_particles_name(std::string name) {
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		particle_name[i] = name; // Default
	}
}

void SystemOfParticles::load_state(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	file >> number_of_particles;
	file >> dt;

	int N = number_of_particles;

	r = dof_mat_t(N, 3);
	v = dof_mat_t(N, 3);
	F = dof_mat_t(N, 3);
	m = Eigen::VectorXd(N);

	particle_name.reset(new std::string[N]);
	f.reset(new Force(true));		//Force of interaction

	for (unsigned int i = 0; i < number_of_particles; i++) {
		file >> particle_name[i] >> m[i] >> r(i, 0) >> r(i, 1) >> r(i, 2) >> v(i, 0) >> v(i, 1) >> v(i, 2);
	}

	file.close();
}

void SystemOfParticles::store_state(std::string file_name) {
	std::ofstream file;
	file.open(file_name);
	file << number_of_particles << std::endl;
	file << dt << std::endl;

	int N = number_of_particles;

	//r = new double[3*N];	//Position
	//v = new double[3*N];	//Velocity
	//F = new double[3*N];	//Force
	//m = new double[N];	//Mass
	//particle_name = new std::string [N];
	//f = new Force();		//Force of interaction

	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		file << particle_name[i] << "\t";
		file << m[i] << "\t";
		file << r(i, 0) << "\t" << r(i, 1) << "\t" << r(i, 2) << "\t";
		file << v(i, 0) << "\t" << v(i, 1) << "\t" << v(i, 2) << std::endl;
	}
	file.close();
}

void SystemOfParticles::store_xyz_file(bool append, std::string file_name, std::string delimiter) {
	std::ofstream file;
	if (append) {
		file.open(file_name, std::ios::app);
	}
	else {
		file.open(file_name);
		file << number_of_particles << std::endl;
	}
	for (unsigned int idx = 0; idx < number_of_particles; idx ++) {
		int i = idx / 3;
		file << particle_name[i];
		file << delimiter << r(idx, 0) << delimiter << r(idx, 1) << delimiter << r(idx, 2) << std::endl;
	}
	file.close();
}


void SystemOfParticles::store_files(int iterator, int factor_xyz, int factor_xy,
	double avg_pressure, double heat_capacity,
	std::string k_energy_file,
	std::string p_energy_file,
	std::string e_energy_file,
	std::string temperature_file,
	std::string average_pressure_file,
	std::string heat_capacity_file) {

	if (iterator == 1) {
		storeData(f->unit_of_energy(), k_energy_file);
		storeData(f->unit_of_energy(), p_energy_file);
		storeData(f->unit_of_energy(), e_energy_file);
		storeData(f->unit_of_temperature(), temperature_file);
		storeData(f->unit_of_pressure(), average_pressure_file);
		storeData(f->unit_of_heat_capacity(), heat_capacity_file);
		store_xyz_file();
	}

	if (!(iterator % factor_xyz)) {
		store_xyz_file(true);
	}
	if (!(iterator % factor_xy)) {

		double time = dt * iterator;
		storeXYData(time, K_energy, k_energy_file);
		storeXYData(time, P_energy, p_energy_file);
		storeXYData(time, energy, e_energy_file);
		storeXYData(time, T, temperature_file);
		storeXYData(time, avg_pressure / time, average_pressure_file);
		storeXYData(time, heat_capacity, heat_capacity_file);
	}
}
void SystemOfParticles::show_infos(int iterator, int factor_ecran, int factor_percent, double average_pressure, double gas_constant) {

	double K, P, E;

	std::cout << std::setprecision(1);

	if (iterator == 1) {
		std::cout << "\n\t========================================\n";
		std::cout << "\t|State\t| K energy  |P energy | Mec. energy | Temperature (K) | Pressure | Gas const.(J/(mol.K)) | Error \t\n\n \t\n\n";
		std::cout.flush();
	}
	if (!(iterator % factor_ecran)) {
		//            if (!((it/fac_percent)%10)) {
		K = knetic_energy();                // alredy computes temperature T
		P = potential_energy();
		E = K + P;
		std::cout << std::fixed << "\t|" << double(iterator / factor_percent) << "%";
		//            std::cout << std::scientific;
		std::cout << "\t|" << K;// << "  " ;
		std::cout << "\t    |" << P;// << "  " ;
		std::cout << "   |" << E;// << "  ";
		std::cout << "\t    |" << T * f->unit_of_temperature();// << "  ";
		std::cout << "\t      |" << average_pressure / (dt * iterator);// << "  ";
		std::cout << "\t |" << gas_constant;
		std::cout << "\t\t\t  |" << abs(gas_constant / f->get_gas_constant() - 1.0) * 100.0 << "%";
		std::cout << std::endl;
		std::cout.flush();
		//            }
	}

}
