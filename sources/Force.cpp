#include "Force.hpp"



Force::Force() : p_F {0, 0, 0} {}

// Lennard-Jonnes Force
// Using Kong rules (CL Kong, The Journal of Chemical Physics. 59 (5): 2464.)
double* Force::operator () (double *r, double *ro, double sigmaA, double sigmaB, double epsilonA, double epsilonB) {
    double Rx, Ry, Rz;
    Rx = r[0] - ro[0];
    Ry = r[1] - ro[1];
    Rz = r[2] - ro[2];
    double Rm2 = 1/(pow(Rx,2.0) + pow(Ry,2.0) + pow(Rz,2.0));
    double Rm6 = pow(Rm2,3);
    double epsilonsigma6A = epsilonA*pow(sigmaA,6.0);
    double epsilonsigma6B = epsilonB*pow(sigmaB,6.0);
    double epsilonsigma6 = sqrt(epsilonsigma6A*epsilonsigma6A);
    double epsilonsigma12 = pow(0.5*(pow(epsilonsigma6A*pow(sigmaA,6),1.0/13.0) + pow(epsilonsigma6B*pow(sigmaB,6),1.0/13.0)),13);
    
    p_unit_space =  (epsilonsigma12/epsilonsigma6)*(epsilonsigma12/epsilonsigma6);
    
    p_unit_energy = (epsilonsigma6*epsilonsigma6/epsilonsigma12);
    p_unit_time = sqrt(atomic_mass/p_unit_energy)*p_unit_space;
    p_unit_temperature = pow(p_unit_space/p_unit_time,2)*atomic_mass;
    
    double C = 24.0*pow(Rm2,4)*(2.0*Rm6 - 1.0);
    
    p_F[0] = C*Rx;///unitForce;
    p_F[1] = C*Ry;///unitForce;
    p_F[2] = C*Rz;///unitForce;
    
    return p_F;
}

double Force::potential(double *r, double *ro, double sigmaA, double sigmaB, double epsilonA, double epsilonB) {
    double Rx, Ry, Rz;
    Rx = r[0] - ro[0];
    Ry = r[1] - ro[1];
    Rz = r[2] - ro[2];
    double Rm2 = 1/(pow(Rx,2.0) + pow(Ry,2.0) + pow(Rz,2.0));
    double Rm6 = pow(Rm2,3);
    
    return 4.0*Rm6*(Rm6 - 1.0);
    
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


