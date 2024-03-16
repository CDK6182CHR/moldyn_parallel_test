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
#include "sycl_helpers.hpp"


#if 0
SystemOfParticles::SystemOfParticles() :
	r(nullptr),
	v(nullptr),
	F(nullptr),
	m(nullptr),
	number_of_particles(0),
	f(nullptr),
	particle_name(nullptr),
	Q(sycl::cpu_selector_v),
	dt(0.0),
	lateral_size(0.0),
	lattice_parameter(0.0)
{}
#endif

SystemOfParticles::SystemOfParticles(int N, double temperature, double dens, double time_step):
	hr(new double[3*N]), hv(new double[3*N]), hF(new double [3*N]), hm(new double[N])
	 ,br(mdutil::as_3d(hr), sycl::range<1>(N)), bv(mdutil::as_3d(hv), sycl::range<1>(N))
	 ,bF(mdutil::as_3d(hF), sycl::range<1>(N)), bm(hm.get(), sycl::range<1>(N))
	, Q(sycl::cpu_selector_v)
 {
	std::cout << "\tSYCL device: " << Q.get_device().get_info<sycl::info::device::name>() << std::endl;
	number_of_particles = N;

	// r.reset(new double[3 * N]);	            //Position
	// v.reset(new double[3 * N]);	            //Velocity
	// F.reset(new double[3 * N]);	            //Force
	// m.reset(new double[N]);	                //Mass
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
	reset_force();
	//for (unsigned int i = 0; i < 3 * number_of_particles; i += 1) {
	//	F[i] = 0.0;
	//}


	std::cout << "\n\tSystem of units used (in S.I. units):\n";
	std::cout << "\tUnit of energy: " << f->unit_of_energy() << std::endl;
	std::cout << "\tUnit of space: " << f->unit_of_space() << std::endl;
	std::cout << "\tUnit of time: " << f->unit_of_time() << std::endl;
	std::cout << "\tUnit of temperature: " << f->unit_of_temperature() << std::endl;

	//Q.wait_and_throw();
	std::cout << "done: set_initial_state" << std::endl;
}

void SystemOfParticles::reset_force()
{
	Q.submit([&, this](sycl::handler& h) {
		auto acc_F = bF.get_access(h, sycl::write_only);
		h.parallel_for(sycl::range(number_of_particles), [=](sycl::id<1> i) {
			std::fill_n(acc_F[i], 3, 0.0);
			});
		});
}

void SystemOfParticles::execute_interations(int number_of_interations) {

	int factor_percent = number_of_interations / 100;
	int factor_ecran = factor_percent;
	int factor_store_state = 1000;
	int factor_xy = 1000;

	double average_temperature = 0.0;
	double average_pressure = 0.0;
	double average_energy = 0.0;
	double square_energy = 0.0;
	double heat_capacity = 0.0;
	double gas_constant = 0.0;

	std::unique_ptr<Timer> timer(new Timer);

	timer->start_timer();
	for (unsigned int it = 1; it < number_of_interations + 1; it += 1) {

		double P_loc = 0, K_loc = 0, U_loc = 0;

#ifdef MANUAL_PARALLEL
#pragma omp parallel reduction(+:P_loc) reduction(+:K_loc) reduction(+:U_loc)
#endif
		{
#ifdef MANUAL_PARALLEL
#pragma omp for
#endif
			reset_force();

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

		//std::cout << "kinetic: " << K_energy << std::endl;

		average_pressure += pressure;
		average_temperature += T;
		average_energy += energy;
		square_energy += energy * energy;
		heat_capacity = (square_energy / it - pow(average_energy / it, 2)) / (T * T);        // From canonical ensemble

		gas_constant = average_pressure * volume / (number_of_particles * average_temperature);
		gas_constant *= f->unit_of_gas_constant();
		
		Q.wait_and_throw();

		store_files(it, factor_store_state, factor_xy, average_pressure, heat_capacity);

		show_infos(it, factor_ecran, factor_percent, average_pressure, gas_constant);

		if (!(it % factor_ecran)) timer->register_time("simulation time: ");

	}
	timer->end_timer();
	//Q.wait_and_throw();
}

void SystemOfParticles::compute_interations() {

	// computes N(N - 1) pairs of interactions 
	// the force array is not set to zero in order to compute
	// the velocity at the last step of the Velocity-Verlet
	// algorithm
#if 0
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for (int i = 0; i < 3 * number_of_particles; i += 3) {
		for (int j = 0; j < 3 * number_of_particles; j += 3) {
			if (i == j) continue;

			auto F_ptr = (*f)(&r[i], &r[j]);   // move construct
			//            for (unsigned int k = 0; k < 3; k += 1) {
			//            
			//                F[i + k] += F_ptr[k];
			//                F[j + k] -= F_ptr[k];
			//                
			//            }

			F[i] += F_ptr[0];
			F[i + 1] += F_ptr[1];
			F[i + 2] += F_ptr[2];

			//F[j] -= F_ptr[0];
			//F[j + 1] -= F_ptr[1];
			//F[j + 2] -= F_ptr[2];

		}
	}
#else

	//sycl::buffer<double[3]> bF(mdutil::as_3d(F), sycl::range<1>(number_of_particles));
	//sycl::buffer<double[3]> br(mdutil::as_3d(r), sycl::range<1>(number_of_particles));
	//try {
		Q.submit([&, this](sycl::handler& h) {
			sycl::accessor aF(bF, h, sycl::read_write);
			sycl::accessor ar(br, h, sycl::read_only);
			h.parallel_for<class ComputeForces>(sycl::range<1>(number_of_particles), [=, this](sycl::id<1> i) {
				for (int j = 0; j < number_of_particles; j++) {
					if (j == i) 
						continue;
					double force_loc[3];
					f->operator()(ar[i], ar[j], force_loc);
					for (int k = 0; k < 3; k++) {
						aF[i][k] += force_loc[k];
					}
				}
				});
			});
		//Q.wait_and_throw();
		//std::cout << "done: compute_interaction" << std::endl;
	//}
	//catch (std::exception& e) {
	//	std::cout << "SYCL error: " << e.what() << std::endl;
	//	std::exit(1);
	//}
#endif
}

void SystemOfParticles::move_particles() {
#if 0
#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < 3 * number_of_particles; i += 3) {
		for (int j = 0; j < 3; j += 1) {

			r[i + j] += v[i + j] * dt + 0.5 * (F[i + j] / m[i / 3]) * dt * dt;

		}
	}
#else
	Q.submit([&, this](sycl::handler& h) {
		auto ar = br.get_access(h, sycl::read_write);
		auto av = bv.get_access(h, sycl::read_only);
		auto aF = bF.get_access(h, sycl::read_only);
		auto am = bm.get_access(h, sycl::read_only);
		h.parallel_for<class MoveParticles>(sycl::range(number_of_particles), [=, this](sycl::id<1> i) {
			for (int k = 0; k < 3; k++) {
				ar[i][k] += av[i][k] * dt + 0.5 * (aF[i][k] / am[i]) * dt * dt;
			}
			});
		});
	//Q.wait_and_throw();
	//std::cout << "done: move_particles" << std::endl;
#endif
}

void SystemOfParticles::compute_velocities() {
#if 0
#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < 3 * number_of_particles; i += 3) {
		for (int j = 0; j < 3; j += 1) {

			v[i + j] += 0.5 * (F[i + j] / m[i / 3]) * dt;

		}
	}
#else 
	Q.submit([&, this](sycl::handler& h) {
		auto av = bv.get_access(h, sycl::read_write);
		auto aF = bF.get_access(h, sycl::read_only);
		auto am = bm.get_access(h, sycl::read_only);
		h.parallel_for(sycl::range(number_of_particles), [=, this](sycl::id<1> i) {
			for (int k = 0; k < 3; k++)
				av[i][k] += 0.5 * (aF[i][k] / am[i]) * dt;
			});
		});
	//Q.wait_and_throw();
	//std::cout << "done: compute_vel" << std::endl;
#endif
}

double SystemOfParticles::check_wall_collisions() {
	// Elastic walls
	double P = 0.0;

#if 0
#ifdef _OPENMP
#pragma omp for 
#endif
	for (int i = 0; i < 3 * number_of_particles; i += 3) {
		for (int j = 0; j < 3; j += 1) {
			if (r[i + j] < 0.0) {
				v[i + j] *= -1.; //- elastic walls
				P += 2.0 * m[i / 3] * abs(v[i + j]) / dt;        // the average forçe from the walls
			}
			if (r[i + j] >= lateral_size) {
				v[i + j] *= -1.;  //- elastic walls
				P += 2.0 * m[i / 3] * abs(v[i + j]) / dt;        // the average forçe from the walls
			}
		}
	}
#else
	// SYCL reduction 
	// spec: https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html#sec:reduction
	// see also: https://github.com/Michoumichmich/MolecularDynamics
	{
		sycl::buffer<double> bP{ &P, sycl::range(1) };
		Q.submit([&, this](sycl::handler& h) {
			auto acc_v = bv.get_access(h, sycl::read_write);
			auto acc_r = br.get_access(h, sycl::read_only);
			auto acc_m = bm.get_access(h, sycl::read_only);
			auto red_P = sycl::reduction(bP, h, sycl::plus<>());
			h.parallel_for(sycl::range(number_of_particles), red_P, [=, this](sycl::id<1> i, auto& p_sum) {
				for (int k = 0; k < 3; k++) {
					if (acc_r[i][k] < 0.0) {
						acc_v[i][k] *= -1.; //- elastic walls
						p_sum += 2.0 * acc_m[i] * std::abs(acc_v[i][k]) / dt;        // the average forçe from the walls
					}
					if (acc_r[i][k] >= lateral_size) {
						acc_v[i][k] *= -1.;  //- elastic walls
						p_sum += 2.0 * acc_m[i] * std::abs(acc_v[i][k]) / dt;        // the average forçe from the walls
					}
				}
				});
			});
	}   // end of reduction region, the destruction of bP should write data into P
#endif

	return P / (6.0 * lateral_size * lateral_size);            // the average forçe form the walls over area

}

double SystemOfParticles::knetic_energy() {
	double v2 = 0.0;

#if 0
#ifdef _OPENMP
#pragma omp for
#endif
	for (int i = 0; i < 3 * number_of_particles; i += 3) {
		for (int j = 0; j < 3; j += 1) {
			v2 += m[i / 3] * v[i + j] * v[i + j];
		}
	}
	//T = v2/(3*(number_of_particles - 1));   // 2024.03.09  for parallel version: put outside
#else 
	{
		sycl::buffer<double> buf_v2{ &v2, 1 };
		Q.submit([&, this](sycl::handler& h) {
			auto acc_m = bm.get_access(h, sycl::read_only);
			auto acc_v = bv.get_access(h, sycl::read_only);
			auto red_v2 = sycl::reduction(buf_v2, h, sycl::plus<>());
			h.parallel_for(sycl::range(number_of_particles), red_v2, [=](sycl::id<1> i, auto& v2_sum) {
				for (int k = 0; k < 3; k++) {
					v2_sum += acc_m[i] * acc_v[i][k] * acc_v[i][k];
				}
				});
			});
	}   // end of reduction
#endif
	return 0.5 * v2;

}

double SystemOfParticles::Temperature() {
	T = 2.0 * knetic_energy() / (3 * number_of_particles);
	return T;
}

double SystemOfParticles::potential_energy() {
	double U = 0.0;
#if 0
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for (int i = 0; i < 3 * number_of_particles - 3; i += 3) {
		for (int j = i + 3; j < 3 * number_of_particles; j += 3) {
			U += f->potential(&r[i], &r[j]);
		}
	}
#else 
	// SYCL reduction
	{
		sycl::buffer<double> buf_U(&U, 1);
		Q.submit([&, this](sycl::handler& h) {
			auto acc_r = br.get_access(h, sycl::read_only);
			auto red_U = sycl::reduction(buf_U, h, sycl::plus<>());
			h.parallel_for(sycl::range(number_of_particles), red_U, [=, this](sycl::id<1> i, auto& sum_U) {
				for (int j = i + 1; j < number_of_particles; j++) {
					sum_U += f->potential(acc_r[i], acc_r[j]);
				}
				});
			});
	}
#endif
	return U;
}

void SystemOfParticles::set_positions() {

	int n;

	n = int(ceil(pow(number_of_particles, 1.0 / 3.0))); // Number of molecules per direction
	lattice_parameter = lateral_size / n;        // Distance between molecules

	
#if 0
	int idx = 0;
	int TN = 3 * number_of_particles;

	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n; j++) {
			for (unsigned int k = 0; k < n; k++) {
				if (idx < TN) {
					r[idx] = (i + 0.5) * lattice_parameter;
					r[idx + 1] = (j + 0.5) * lattice_parameter;
					r[idx + 2] = (k + 0.5) * lattice_parameter;
				}
				idx += 3;
			}
		}
	}
#else
	// For simplicity, just use single_task ...
	Q.submit([&, this](sycl::handler& h) {
		auto acc_r = br.get_access(h, sycl::write_only);
		h.single_task([=, this]() {
			int idx = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						if (idx < number_of_particles) {
							acc_r[idx][0] = (i + 0.5) * lattice_parameter;
							acc_r[idx][1] = (j + 0.5) * lattice_parameter;
							acc_r[idx][2] = (k + 0.5) * lattice_parameter;
						}
						idx++;
					}
				}
			}
			});
		});
#endif
	//Q.wait_and_throw();
	//std::cout << "done: set_pos" << std::endl;
}

void SystemOfParticles::set_velocities(double mean, double dispertion) {

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	Q.submit([&, this](sycl::handler& h) {
		auto acc_v = bv.get_access(h, sycl::read_write);
		auto acc_m = bm.get_access(h, sycl::read_only);
		h.single_task([=, this]() {
			double Vcm[3] = { 0.,0.,0. };
			double M = 0.0;
			double v2 = 0.0;
			// construct a trivial random generator engine from a time-based seed:
			
			std::default_random_engine generator(seed);
			std::normal_distribution<double> distribution(mean, dispertion);

			for (int i = 0; i < number_of_particles; i++) {
				for (int k = 0; k < 3; k++)
					acc_v[i][k] = distribution(generator);
			}
			for (int i = 0; i < number_of_particles; i++) {
				for (int j = 0; j < 3; j += 1) {
					Vcm[j] += acc_m[i] * acc_v[i][j];
				}
			}

			for (unsigned int i = 0; i < number_of_particles; i += 1) {
				M += acc_m[i];
			}

			for (unsigned int j = 0; j < 3; j += 1) {
				Vcm[j] /= number_of_particles * M;
			}

			for (unsigned int i = 0; i < number_of_particles; i++) {
				for (unsigned int j = 0; j < 3; j += 1) {
					acc_v[i][j] -= Vcm[j];
					v2 += acc_m[i] * acc_v[i][j] * acc_v[i][j];
				}
			}

			double To = v2 / (3 * (number_of_particles - 1));
			double temperature_factor = sqrt(T / To);

			for (int i = 0; i < number_of_particles; i++) {
				for (int j = 0; j < 3; j += 1) {
					acc_v[i][j] *= temperature_factor;
				}
			}
			});
		});
	//Q.wait_and_throw();
	//std::cout << "done: set_vel" << std::endl;
}

void SystemOfParticles::set_equal_masses(double mass) {
#if 0
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		m[i] = mass;
	}
#else 
	Q.submit([&](sycl::handler& h) {
		auto acc_m = bm.get_access(h, sycl::write_only);
		h.parallel_for(sycl::range(number_of_particles), [=](sycl::id<1> i) {
			acc_m[i] = mass;
			});
		});
	//Q.wait_and_throw();
	std::cout << "done: set_equal_mass" << std::endl;
#endif
}

void SystemOfParticles::set_particles_name(std::string name) {
	for (unsigned int i = 0; i < number_of_particles; i += 1) {
		particle_name[i] = name; // Default
	}
}

#if 0
void SystemOfParticles::load_state(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	file >> number_of_particles;
	file >> dt;

	int N = number_of_particles;

	r.reset(new double[3 * N]);	//Position
	v.reset(new double[3 * N]);	//Velocity
	F.reset(new double[3 * N]);	//Force
	m.reset(new double[N]);	//Mass
	particle_name.reset(new std::string[N]);
	f.reset(new Force(true));		//Force of interaction

	for (unsigned int i = 0; i < 3 * number_of_particles; i += 3) {
		file >> particle_name[i / 3] >> m[i / 3] >> r[i] >> r[i + 1] >> r[i + 2] >> v[i] >> v[i + 1] >> v[i + 2];
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
		file << r[i] << "\t" << r[i + 1] << "\t" << r[i + 2] << "\t";
		file << v[i] << "\t" << v[i + 1] << "\t" << v[i + 2] << std::endl;
	}
	file.close();
}
#endif

void SystemOfParticles::store_xyz_file(bool append, std::string file_name, std::string delimiter) {
	std::ofstream file;
	if (append) {
		file.open(file_name, std::ios::app);
	}
	else {
		file.open(file_name);
		file << number_of_particles << std::endl;
	}
	auto acc_r = br.get_host_access(sycl::read_only);

	for (int idx = 0; idx < number_of_particles; idx++) {
		file << particle_name[idx];
		file << delimiter << acc_r[idx][0] << delimiter << acc_r[idx][1] << delimiter << acc_r[idx][2] << std::endl;
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
