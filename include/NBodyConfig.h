
/**
 * A header to control the configuration and compilation
 * of all the N-Body code.
 */

#ifndef _NBODY_CONFIG_
#define _NBODY_CONFIG_


/**
 * Init the points to be two Plummer clusters.
 */
static const int NBODY_INIT_TWO_CLUSTERS = 1;

/**
 * Distance between center of two clusters in each dimension.
 */
static const double NBODY_TWO_CLUSTER_SEP = 2.0;

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
static const int NBODY_NPROCS = 2;
#else
static const int NBODY_NPROCS = 1;
#endif

static const int SORT_EXCHANGE_SIZE = 64;

static const int DEFAULT_HASHED_OCTREE_H = 14;


#endif
