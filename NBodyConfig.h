
/**
 * A header to control the configuration and compilation
 * of all the N-Body code.
 */

#ifndef _NBODY_CONFIG_
#define _NBODY_CONFIG_


/**
 * Should the renderer be built and run with the simulation?
 */
#define NBODY_SIM_WITH_RENDERER 1

/**
 * Should the simulation be built to run in parallel?
 */
#define NBODY_PARALLEL 1

/**
 * Init the points to be two Plummer clusters.
 */
static const int NBODY_INIT_TWO_CLUSTERS = 1;

/**
 * Init the points to be one Plummer clusters.
 */
static const int NBODY_INIT_PLUMMER = 0;

/**
 * Init the points to be uniformly distributed in unit sphere.
 */
static const int NBODY_INIT_UNIFORM = 0;

/**
 * The number of processors the simulation should use.
 */
#if NBODY_PARALLEL
static const int NBODY_NPROCS = 8;
#else
static const int NBODY_NPROCS = 1;
#endif

#endif
