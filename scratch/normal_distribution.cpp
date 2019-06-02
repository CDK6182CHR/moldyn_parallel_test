// normal_distribution
#include <iostream>
#include <string>
#include <chrono>
#include <random>

double gaussdist();
 
int main()
{
//   construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);

  std::normal_distribution<double> distribution (0.0,40.0);

  std::cout << "some Normal-distributed(0.0,1.0) results:" << std::endl;
  for (int i=0; i<10; ++i)
    std::cout << distribution(generator) << std::endl;

  std::cout << "gaussdist() results:" << std::endl;
  for (int i=0; i<10; ++i)
    std::cout << gaussdist() << std::endl;

  return 0;
// Bad Idea

}

//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
  static bool available = false;
  static double gset;
  double fac, rsq, v1, v2;
  if (!available) {
  do {
    v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
    v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
    rsq = v1 * v1 + v2 * v2;
  } while (rsq >= 1.0 || rsq == 0.0);

  fac = sqrt(-2.0 * log(rsq) / rsq);
  gset = v1 * fac;
  available = true;

  return v2*fac;
  } else {

  available = false;
  return gset;

}
}

//int main()
//{
//  const int nrolls=10000;  // number of experiments
//  const int nstars=100;    // maximum number of stars to distribute

//  std::default_random_engine generator;
//  std::normal_distribution<double> distribution(7.0,2.0);

//  int p[10]={};

//  for (int i=0; i<nrolls; ++i) {
//    double number = distribution(generator);
//    if ((number>=0.0)&&(number<10.0)) ++p[int(number)];
//  }

//  std::cout << "normal_distribution (5.0,2.0):" << std::endl;

//  for (int i=0; i<10; ++i) {
//    std::cout << i << "-" << (i+1) << ": ";
//    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
//  }

//  return 0;
//}
