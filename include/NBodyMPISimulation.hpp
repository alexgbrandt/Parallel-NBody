
#ifndef _NBODY_MPI_SIM_
#define _NBODY_MPI_SIM_

/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation is serial and performs the exact, but naive N^2 algorithm
 * for force calculation.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 */
void NBodyNaive_1(long N, double dt, double t_end, time_t seed);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the naive N^2 algorithm
 * but evenly distributes the bodies across MPI ranks
 * for parallel force claculation and position update.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 */
void NBodyNaive_NBMPI(long N, double dt, double t_end, time_t seed);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step the root process builds the octree,
 * broadcasts it, and then each processor computes
 * but evenly distributes the bodies across MPI ranks
 * for parallel force claculation and position update.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH1_NBMPI(long N, double dt, double t_end, time_t seed, double theta);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Trees are iteratively merged
 * in a reduce, with the final global tree then broadcasted.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH2_NBMPI(long N, double dt, double t_end, time_t seed, double theta);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Trees are iteratively merged
 * in a reduce, with the final global tree then broadcasted.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * Here, the bodies are sorted at each step on the root and then distributed
 * across ranks using their Morton ordering to achieve a spatial decomposition.
 * This makes the octree merges nearly instant.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH3_NBMPI(long N, double dt, double t_end, time_t seed, double theta);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Trees are iteratively merged
 * in a reduce, with the final global tree then broadcasted.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * Here, the bodies are sorted at each step on the root and then distributed
 * accross ranks using their Morton ordering to achieve a spatial decomposition.
 * Further, bodies are distributed unevenly, but based on a work-estimate from
 * the previous simulation step in order to achieve load-balance.
 *
 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH4_NBMPI(long N, double dt, double t_end, time_t seed, double theta);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Trees are iteratively merged
 * in a reduce, with the final global tree then broadcasted.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * Here, the bodies are sorted at each step using a distributed sort
 * using their Morton ordering to achieve a spatial decomposition.
 * After the sort, bodies are exchanged between ranks based on a work-estimate from
 * the previous simulation step in order to achieve load-balance.
 * This avoids the gather to root of NBodyBH4 for sorting.

 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH5_NBMPI(long N, double dt, double t_end, time_t seed, double theta);


/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Any data needed from non-local parts
 * of the global octree is dynamically requested during the tree traversal.
 * This avoids the requirement to explictly share the global tree
 * at every step with every processor.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * Here, the bodies are sorted at each step using a distributed sort
 * using their Morton ordering to achieve a spatial decomposition.
 * After the sort, bodies are exchanged between ranks based on a work-estimate from
 * the previous simulation step in order to achieve load-balance.
 * This avoids the gather to root of NBodyBH4 for sorting.

 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH6_NBMPI(long N, double dt, double t_end, time_t seed, double theta);



//distributed spatial decomp and load balancing
//build local tree, broadcast branch nodes,
//dynamic and asynchronous request data during traversal
/**
 * Run a N-body gravity simulation with t_end/dt time steps,
 * and a random seed.
 *
 * This simulation uses the Barnes-Hut NlogN octree algorithm.
 * At every step each rank builds a "local octree"
 * from its local bodies. Any data needed from non-local parts
 * of the global octree is dynamically and asynchronously requested
 * during the tree traversal. Asynchronous communication allows
 * the tree traversal to continue without waiting for the
 * response for non-local data.
 * Each rank then computes forces and position updates for its local bodies.
 *
 * Here, the bodies are sorted at each step using a distributed sort
 * using their Morton ordering to achieve a spatial decomposition.
 * After the sort, bodies are exchanged between ranks based on a work-estimate from
 * the previous simulation step in order to achieve load-balance.
 * This avoids the gather to root of NBodyBH4 for sorting.

 * @param N: number of bodies in the simulation
 * @param dt: the timestep
 * @param t_end: the total simulation time
 * @param seed: the random seed used for initialization
 * @param theta: paramreter controller MAC for the B-H force calculation.
 */
void NBodyBH7_NBMPI(long N, double dt, double t_end, time_t seed, double theta);

#endif