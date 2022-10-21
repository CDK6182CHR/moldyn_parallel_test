# Moldyn

A particle dynamics simulator

## GENERAL NOTES

- This is a 3D particle simulation program for general purporse (`main.cpp` contains a sample usage).
- Tested in Linux (Ubuntu 18.04, 20.04) 64bits system.
- Clone or download, open terminal on folder, `make` and then `make run`.
- You can set parameter in the main.cpp file.
=======

## BASIC GUIED

- Moldyn program simulates a 3D system of interacting particles via [Lennard-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential). The default simulation system is made of Argon atoms enclosed in a elastic wall recipient (refletions on the wall are perfectly elastic), therefore the mechanical energy is constant ([microcanonical ensemble](https://en.wikipedia.org/wiki/Microcanonical_ensemble)), while the temperature fluctuates around an average value, and the system undergoes throught sussesives isochoric themodynamical transformations untill equilibrium as the simulation converges. The convergence can be live checked as the "Error" outuput on screen, which computs the relative error on measuring the [gas constant](https://en.wikipedia.org/wiki/Gas_constant).

- The inputs are (check `main.cpp` file):
    1. Number of particles (works better with powers of 3).
    2. The density of the system in J.mol-1 units.
    3. The initial temperature in K.
    
- The default data is for a system of Argon atoms. Run the default simulation once to check on the screen the system of units used.

- The outuput files are:
    1. `Avg_Pressure.dat`, which stores the average pressure of the system with time.
    2. `Enegy.dat`, containing the total energy as a function of time, and the same for knetic and potential energies, respectively:
    3. `kEnergy.dat`.
    4. `pEnergy.dat`.
    5. `heat_capacity.dat`: althought the system is better classified as a microcanonical ensamble, the [isochoric heat capacity](https://en.wikipedia.org/wiki/Isochoric_process) is computed since, due to numerical convergence errors, the system's energy fluctuates and can be treated as a [canonical ensamble](https://en.wikipedia.org/wiki/Canonical_ensemble).
    6. `positions.xyz` is a file in [xyz format](https://en.wikipedia.org/wiki/XYZ_file_format), and can be read with [VMD](https://www.ks.uiuc.edu/Research/vmd/).
    7. `Temperature.dat` contains the temperature as a function of time.
=======
