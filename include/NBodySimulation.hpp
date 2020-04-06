

#ifndef _NBODY_SIMULTATION_
#define _NBODY_SIMULTATION_

/**
 * Start the gravitational N-Body simultation in serial.
 *
 * @param N, the number of bodies.
 * @param dt, the time step.
 * @param t_end, the total time to simulation.
 * @param seed, a random number seed, 0 to use the current time as seed.
 * @param theta, a parameter for the Barnes-Hut multipole acceptance criteria.
 */
void NBodySimSerial(long N, double dt, double t_end, time_t seed, double theta);


#if NBODY_PARALLEL
/**
 * Start the gravitational N-Body simultation in parallel.
 *
 * @param N, the number of bodies.
 * @param dt, the time step.
 * @param t_end, the total time to simulation.
 * @param seed, a random number seed, 0 to use the current time as seed.
 * @param theta, a parameter for the Barnes-Hut multipole acceptance criteria.
 */
void NBodySimParallel(long N, double dt, double t_end, time_t seed, double theta);
#endif

#endif