#include "Force.hpp"
#include <iostream>



Force::Force(double sigma, double epsilon, double mass)
{

    p_unit_energy *= epsilon;                                                      // In units of 993.65322254 J.mol-1 (for Ar)
    p_unit_space *= sigma;                                                         // In units of Borh Radius
    p_unit_mass *= mass;                                                           // In units of atomic mass
    p_unit_time = sqrt(p_unit_mass/p_unit_energy)*p_unit_space;
    p_unit_temperature = pow(p_unit_space,2)*p_unit_mass/pow(p_unit_time,2)/Kb;
    p_unit_pressure = p_unit_energy/pow(p_unit_space,3);
    p_unit_gas_constant = (p_unit_energy/p_unit_temperature)*avogadro_constant;       // J.mol-1 
    p_unit_heat_capacity = pow(p_unit_energy,2)/pow(p_unit_temperature,2); 
}


vec3_t Force::operator () (const Eigen::Ref<vec3_t>& r, const Eigen::Ref<vec3_t>& ro)const {
#if 0
    vec3_t force{};

    double Rx = r[0] - ro[0];
    double Ry = r[1] - ro[1];
    double Rz = r[2] - ro[2];
    
    double Rm2 = 1.0/(Rx*Rx + Ry*Ry + Rz*Rz);
    double Rm6 = Rm2*Rm2*Rm2;
    
    double C = 48.0*pow(Rm2,4)*(Rm6 - 0.5);

    force[0] = C*Rx;///unitForce;
    force[1] = C*Ry;///unitForce;
    force[2] = C*Rz;///unitForce;
    
    return force;
#else 
    auto dr = r - ro;
    double Rm2 = 1.0 / (dr.matrix().squaredNorm());
    double Rm6 = std::pow(Rm2, 3);
    double C = 48.0 * std::pow(Rm2, 4) * (Rm6 - 0.5);
    return C * dr;
#endif
}

double Force::potential(const Eigen::Ref<vec3_t>& r, const Eigen::Ref<vec3_t>& ro) {
#if 0
    double Rx = r[0] - ro[0];
    double Ry = r[1] - ro[1];
    double Rz = r[2] - ro[2];

    double Rm2 = 1.0/(Rx*Rx + Ry*Ry + Rz*Rz);
    double Rm6 = Rm2*Rm2*Rm2;

    return 4.0*Rm6*(Rm6 - 1.0);
#else
    auto dr = r - ro;
    double Rm2 = 1.0 / (dr.matrix().squaredNorm());
    double Rm6 = std::pow(Rm2, 3);
    return 4.0 * Rm6 * (Rm6 - 1.0);
#endif
}

double Force::unit_of_time() {
    return p_unit_time;
}
double Force::unit_of_space() {
    return p_unit_space;
}
double Force::unit_of_energy() {
    return p_unit_energy;
}
double Force::unit_of_temperature() {
	return p_unit_temperature;
}
double Force::unit_of_pressure(){
    return p_unit_pressure;
}
double Force::unit_of_heat_capacity() {
    return p_unit_heat_capacity;
}
double Force::unit_of_gas_constant() {
    return p_unit_gas_constant;
}
double Force::get_avogadro_constant() {
    return avogadro_constant;
}
double Force::get_boltzman_constant() {
    return Kb;
}
double Force::get_gas_constant() {
    return gas_constant;
}
