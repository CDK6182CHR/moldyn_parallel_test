#ifndef FORCELJ_H
#define FORCELJ_H
#include <math.h>

// Lennard-Jonnes Potential for modeling interactiong particles

class Force
{
    private:
        
        constexpr static double Kb = 1.38064852e-23;                                           // Boltzmann constant in m²kg/s²K.
        constexpr static double atomic_mass = 1.66053906660e-27;                               // 1/12 of C12's mass in kg.
        constexpr static double bohr_radius = 5.2917721067e-11;                                // in m.
        constexpr static double elementary_charge = 1.6e-19;                                   // 1eV for unit of energy in J.
        constexpr static double avogadro_constant = 6.02214076e23;                             // for unit of energy in J.mol-1.
        constexpr static double gas_constant = avogadro_constant*Kb;
        constexpr static double sigma_Ar = 3.4e-10/bohr_radius;                                // For Argon in bohr_radius units
        constexpr static double epsilon_Ar = 1.65e-21/(993.65322254/avogadro_constant);        // For Argon in 993.65 J.mol-1 units
        constexpr static double mass_Ar = 39.948;                                              // For Argon in units of atomic mass
                
        double p_unit_mass = mass_Ar*atomic_mass;                                              // Atomic mass of Ar
        double p_unit_energy = 993.65322254/avogadro_constant;                                 // ~ 993.65 J.mol-1 for Ar
        double p_unit_space = bohr_radius;
        double p_unit_time = sqrt(p_unit_mass/p_unit_energy)*p_unit_space;                     // = 5.390954795464439e-15s 
        double p_unit_temperature = pow(p_unit_space,2)*p_unit_mass/pow(p_unit_time,2)/Kb;     // times Kb becomes K
        double p_unit_pressure = p_unit_energy/pow(p_unit_space,3);
        double p_unit_gas_constant = (p_unit_energy/p_unit_temperature)*avogadro_constant;       // J.mol-1
        double p_unit_heat_capacity = pow(p_unit_energy,2)/pow(p_unit_temperature,2);
        
        double p_force[3];
        
    public:
        Force(double sigma = sigma_Ar, double epsilon = epsilon_Ar, double mass = mass_Ar);

        double* operator () (double *r, double *ro);
                
        double potential(double *r, double *ro);
        
        double unit_of_time();
        double unit_of_energy();
        double unit_of_temperature();
        double unit_of_space();
        double unit_of_pressure();
        double unit_of_heat_capacity();
        double unit_of_gas_constant();
        double get_avogadro_constant();
        double get_boltzman_constant();
        double get_gas_constant();
//        std::string unit_of_time_name();
//        std::string unit_of_energy_name();
//        std::string unit_of_temperature_name();
//        std::string unit_of_space_name();
        
};
#endif
