#ifndef FORCELJ_H
#define FORCELJ_H
// Lennard-Jonnes Potential
#include <math.h>
//#include <string>

class Force
{
    private:
        double p_F[3];
        constexpr static double Kb = 1.38064852e-23;                                           // Boltzmann constant in m²kg/s²K.
        constexpr static double atomic_mass = 1.66053906660e-27;                               // 1/12 of C12's mass in kg.
        constexpr static double bohr_radius = 5.2917721067e-11;                                // in m.
        constexpr static double elementary_charge = 1.6e-19;                                   // 1eV for unit of energy in J.
        constexpr static double avogadro_constant = 6.02214076e23;                             // for unit of energy in J.mol-1.
        double p_unit_mass = 39.948*atomic_mass;                                               // Atomic mass of Ar
        double p_unit_energy = 993.65322254/avogadro_constant;                                 // J.mol-1 for Ar
        double p_unit_space = bohr_radius;
        double p_unit_time = sqrt(p_unit_mass/p_unit_energy)*p_unit_space;                     // = 5.390954795464439e-15 
        double p_unit_temperature = pow(p_unit_space,2)*p_unit_mass/pow(p_unit_time,2)/Kb;     // times Kb becomes K
        constexpr static double sigma_Ar = 3.4e-10/bohr_radius;                                // For Argon in bohr_radius units
        constexpr static double epsilon_Ar = 1.65e-21/elementary_charge;                       // For Argon in elementary_charge units
        
        
    public:
        Force();
        
        double* operator () (double *r, double *ro,
                            double sigmaA = sigma_Ar, double sigmaB = sigma_Ar,
                            double epsilonA = epsilon_Ar, double epsilonB = epsilon_Ar);
        
        double potential(double *r, double *ro,
                         double sigmaA = sigma_Ar, double sigmaB = sigma_Ar,
                         double epsilonA = epsilon_Ar, double epsilonB = epsilon_Ar);
        
        double unit_of_time();
        double unit_of_energy();
        double unit_of_temperature();
        double unit_of_space();

//        std::string unit_of_time_name();
//        std::string unit_of_energy_name();
//        std::string unit_of_temperature_name();
//        std::string unit_of_space_name();
        
};
#endif
