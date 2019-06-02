#ifndef SYSTEMOFPARTICLES_H
#define SYSTEMOFPARTICLES_H

class SystemOfParticles
{

    private:
        double *r, *v, *m;
    
    public:
    	SystemOfParticles ();
        SystemOfParticles (int N);
};

#endif
