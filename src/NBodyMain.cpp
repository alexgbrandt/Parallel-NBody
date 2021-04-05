
#include <stdlib.h>

#include "NBodyConfig.h"
#include "NBodySimulation.hpp"

#if NBODY_MPI
#include "NBodyMPISimulation.hpp"
#endif


#include "NBodySerialize.h"
#include <stdlib.h>

/**
 * Main method to get command line args and start and the NBody sim.
 */
int main(int argc, char** argv) {
    long N = 10;
    double dt = 0.01;
    double t_end = 10.0;
    time_t seed = 0;
    double theta = 0.5;
    int algChoice = 7;

    if (argc > 1 && atol(argv[1]) > 0) {
        N = atol(argv[1]);
    }
    if (argc > 2 && atof(argv[2]) > 0) {
        dt = atof(argv[2]);
    }
    if (argc > 3 && atof(argv[3]) > 0) {
        t_end = atof(argv[3]);
    }
    if (argc > 4 && atof(argv[4]) > 0) {
        theta = atof(argv[4]);
    }
    if (argc > 5 && atol(argv[5]) > 0) {
        seed = atol(argv[5]);
    }
    if (argc > 6 && atol(argv[6]) >= 0) {
        algChoice = atol(argv[6]);
    }

#if NBODY_MPI
    switch(algChoice) {
        case 0 :
            NBodyNaive_NBMPI(N, dt, t_end, seed);
            break;
        case 1 :
            NBodyBH1_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 2 :
            NBodyBH2_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 3 :
            NBodyBH3_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 4 :
            NBodyBH4_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 5 :
            NBodyBH5_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 6 :
            NBodyBH6_NBMPI(N, dt, t_end, seed, theta);
            break;
        case 7 :
            NBodyBH7_NBMPI(N, dt, t_end, seed, theta);
            break;
        default :
            break;
    }

#elif NBODY_PARALLEL
    NBodySimParallel(N, dt, t_end, seed, theta);
#else
    NBodySimSerial(N, dt, t_end, seed, theta);
#endif

    return 0;
}