
#include <stdlib.h>

#include "NBodyConfig.h"
#include "NBodySimulation.hpp"

/**
 * Main method to get command line args and start and the NBody sim.
 */
int main(int argc, char** argv) {
    long N = 10;
    double dt = 0.01;
    double t_end = 10.0;
    time_t seed = 0;
    double theta = 0.5;

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

#if NBODY_PARALLEL
    NBodySimParallel(N, dt, t_end, seed, theta);
#else
    NBodySimSerial(N, dt, t_end, seed, theta);
#endif
}