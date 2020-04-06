
#ifndef _NBODY_PARALLEL_
#define _NBODY_PARALLEL_


#include "parallel/include/ExecutorThreadPool.hpp"

#include "NBodyConfig.h"
#include "NBodyOctree.h"

#include <functional>


/**
 * Given the position, velocity, work, (and color) arrays,
 * sort them based on the index map which resulted from sorting the keys.
 *
 * @param N, the number of bodies
 * @param r, an array of 3N values for position.
 * @param v, an array of 3N values for velocitites.
 * @param work, an array of N values for work estimates.
 * @param colors, an array of 4N values for rendering colors of the bodies.
 * @param idx, the index map such that idx[i] should move to i.
 * @param tmpList, a few temporary working space lists, each of size N.
 */
void parallelSortByIdxMap_NB(
	long N,
	double* r,
	double* v,
	double* work,
#if NBODY_SIM_WITH_RENDERER
	float* colors,
#endif
	long* idx,
	long* tmpList[NBODY_NPROCS]);

/**
 * Using the map-reduce pattern, build an octree in parallel using
 * NBODY_NPROCS extra processors. The arrays startN and numN describes the
 * subdata each thread should be responsible for.
 *
 * @param r, an array of 3N values for position.
 * @param m, an array of N values for mass.
 * @param domainSize, the size of the entire octree.
 * @param[in,out] trees, an array of octrees which can be re-used for this process.
 *                       the final tree is also returned in trees[0].
 * @param startN, an array of the starting points for each data partition.
 * @param numN, the size of each data partition.
 *
 * @return 1 iff the tree was successfully built inplace (w.r.t trees[0]).
 */
int mapReduceBuildOctreesInPlace_NB(
	const double* __restrict__ r,
	const double* __restrict__ m,
	double domainSize,
	NBOctree_t* trees[NBODY_NPROCS],
	long startN[NBODY_NPROCS],
	long numN[NBODY_NPROCS]);


/**
 * Compute the acceleration for n bodies at positions r with mass m
 * making use of the octree and a Barnes-Hut MAC for multipole approximation.
 * If a node is internal, the quadrupole moments are included in the potential
 * calculation.
 * This function executed in parallel use data partitions specified
 * by startN and numN.
 *
 * @param n, the number of bodies
 * @param m, an array of n doubles holding the masses of the bodies
 * @param r, an array of 3*n doubles holding the positions of the bodies
 * @param[out] a, an array of 3*n doubles to hold the resulting acceleration of the bodies.
 * @param tree, an octree holding the n bodies.
 * @param list1, an array used as working space of size at least N.
 * @param list2, an array used as working space of size at least N.
 * @param thetaMac, a MAC parameter to control the use of direct calculation or multipole approximation.
 * @param startN, an array of the starting points for each data partition.
 * @param numN, the size of each data partition.
 */
void computeForcesOctreeBHParallel_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* work,
	double* __restrict__ a,
	const NBOctree_t* tree,
	const NBOctreeNode_t** list1[NBODY_NPROCS],
	const NBOctreeNode_t** list2[NBODY_NPROCS],
	double thetaMAC,
	long startN[NBODY_NPROCS],
	long numN[NBODY_NPROCS]);


/**
 * Perform the first half of the leafprog integration in,
 * "kick, drift", in parallel. Updating velocities to the half step
 * and positions to the full step.
 *
 * @param n, the number of bodies.
 * @param dt, the time step for integration.
 * @param r, an array of 3*n doubles holding the positions of the bodies.
 * @param v, an array of 3*n doubles holding the velocitites of the bodies.
 * @param a, an array of 3*n doubles holding the acceleration of the bodies.
 * @param m, an array of n doubles holding the masses of the bodies.
 * @param startN, an array of the starting points for each data partition.
 * @param numN, the size of each data partition.
 */
void performNBodyHalfStepAParallel_NB(
	double dt,
	double* __restrict__ r,
	double* __restrict__ v,
	const double* __restrict__ a,
	const double* __restrict__ m,
	long startN[NBODY_NPROCS],
	long numN[NBODY_NPROCS]);


/**
 * Perform the second half of the leafprog integration,
 * "kick2", in parllel. Updating velocities to the full step
 * from the half step
 *
 * @param n, the number of bodies.
 * @param dt, the time step for integration.
 * @param r, an array of 3*n doubles holding the positions of the bodies.
 * @param v, an array of 3*n doubles holding the velocitites of the bodies.
 * @param a, an array of 3*n doubles holding the acceleration of the bodies.
 * @param m, an array of n doubles holding the masses of the bodies.
 * @param startN, an array of the starting points for each data partition.
 * @param numN, the size of each data partition.
 */
void performNBodyHalfStepBParallel_NB(
	double dt,
	const double* __restrict__ r,
	double* __restrict__ v,
	const double* __restrict__ a,
	const double* __restrict__ m,
	long startN[NBODY_NPROCS],
	long numN[NBODY_NPROCS]);





#endif