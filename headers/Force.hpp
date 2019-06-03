#ifndef FORCELJ_H
#define FORCELJ_H
// Lennard-Jonnes Potential
#include <math.h>
class Force
{
    private:
        double p_F[3];
        double const Kb = 1.38064852e-23; // Boltzmann constant in m²kg/s²K
        double const atomic_mass = 1.66053906660e-27; // 1/12 of C12's mass in kg
        double const bohr_radius = 5.2917721067e-11; // in m
        double p_unit_time;// = bohr_radius/speed_of_light;
        double p_unit_space;
        double p_unit_energy;// = bohr_radius*pow(atomic_mass,2)/pow(p_unit_time,2);
        double p_unit_temperature;// = Kb*pow(p_unit_time,2)/pow(bohr_radius,2)/atomic_mass;
        
        
    public:
        Force();
        
        double* operator () (double *r, double *ro, double RminA = 1.0, double RminB = 1.0, double epsilonA = 1.0, double epsilonB = 1.0);
        
        double potential(double *r, double *ro, double RminA = 1.0, double RminB = 1.0, double epsilonA = 1.0, double epsilonB = 1.0);
        
        double unit_of_time();
        double unit_of_energy();
        double unit_of_temperature();
        double unit_of_space();
        
};
#endif
