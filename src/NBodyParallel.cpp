

#include "parallel/include/ExecutorThreadPool.hpp"

#include "NBodyConfig.h"
#include "NBodyOctree.h"
#include "NBodyForces.h"
#include "NBodyIntegrator.h"
#include "NBodyKeys.h"

#include <functional>

#include <string.h>

/**
 * Implements the reduce parallel programmign pattern using a thread pool
 * and an array of trees to merge. Merging is done pairwise.
 */
void _reduceOctrees_NB(ExecutorThreadPool& threadPool, NBOctree_t** trees, int ntrees) {

	//get steps by bit-hacked ceil(log_2(ntrees))
	int nsteps = 0;
	int ntreetmp = ntrees;
	while (ntreetmp >>= 1) { ++nsteps; }
	if ((ntrees & (ntrees-1))) { ++nsteps; } //round up if not exactly a power of 2

	int stepSize = 1;
	for (int i = 0; i < nsteps; ++i) {
		for (int k = 0; k < ntrees; k += 2*stepSize) {
			if (k + stepSize < ntrees) {
				// fprintf(stderr, "(%d, %d) merge\n", k, k+stepSize);
				std::function<void()> f = std::bind(mergeOctreesInPlace_NB, trees[k], trees[k+stepSize]);
				threadPool.addTaskAtIdx(f, k);
			}
		}
		threadPool.waitForAllThreads(); //sync before stepping;
		stepSize <<= 1;
	}
}


/**
 * Simple function to wrap the build of an Octree to make it void.
 * @see buildOctreeInPlace_NB
 */
void _buildOctreeInPlaceVoid_NB(long N, const double* r, const double* m, double domainSize, NBOctree_t** tree_ptr, int* ret) {
	*ret = buildOctreeInPlace_NB(N, r, m, domainSize, tree_ptr);
	return;
}


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
	long* tmpList[NBODY_NPROCS])
{

#if NBODY_SIM_WITH_RENDERER
	if (NBODY_NPROCS >= 4) {
#else
	if (NBODY_NPROCS >= 3) {
#endif
		ExecutorThreadPool& threadPool = ExecutorThreadPool::getThreadPool();
		std::function<void()> f0 = [=]() {
			memcpy(tmpList[0], idx, sizeof(long)*N);
        	sortByIdxMap3N_NB(N, tmpList[0], r);
		};

		std::function<void()> f1 = [=]() {
			memcpy(tmpList[1], idx, sizeof(long)*N);
        	sortByIdxMap3N_NB(N, tmpList[1], v);
		};
		std::function<void()> f2 = [=]() {
			memcpy(tmpList[2], idx, sizeof(long)*N);
        	sortByIdxMap_NB(N, tmpList[2], work);
		};
		threadPool.addTaskAtIdx(f0, 0);
		threadPool.addTaskAtIdx(f1, 1);
		threadPool.addTaskAtIdx(f2, 2);

#if NBODY_SIM_WITH_RENDERER
		std::function<void()> f3 = [=]() {
			memcpy(tmpList[3], idx, sizeof(long)*N);
        	sortByIdxMap4N_NB(N, tmpList[3], colors);
		};
		threadPool.addTaskAtIdx(f3, 3);
#endif

		threadPool.waitForAllThreads();

	} else {
        memcpy(tmpList[0], idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, tmpList[0], r);
        memcpy(tmpList[0], idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, tmpList[0], v);
        memcpy(tmpList[0], idx, sizeof(long)*N);
        sortByIdxMap_NB(N, tmpList[0], work);
#if NBODY_SIM_WITH_RENDERER
        memcpy(tmpList[0], idx, sizeof(long)*N);
        sortByIdxMap4N_NB(N, tmpList[0], colors);
#endif

	}
}


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
	NBOctree_t** trees,
	long startN[NBODY_NPROCS],
	long numN[NBODY_NPROCS])
{
	ExecutorThreadPool& threadPool = ExecutorThreadPool::getThreadPool();

	int retVals[NBODY_NPROCS];

	//map
	for (int i = 0; i < NBODY_NPROCS; ++i) {
		std::function<void()> f = std::bind(_buildOctreeInPlaceVoid_NB,
			numN[i],
			r + 3*startN[i],
			m + startN[i],
			domainSize,
			trees + i,
			retVals + i
		);
		threadPool.addTaskAtIdx(f, i);
	}

	threadPool.waitForAllThreads(); //sync


	fprintf(stderr, "OCTREES\n\n");
	fprintf(stderr, "Octree[0]\n");
	printOctree_NB(trees[0]);
	fprintf(stderr, "---------------\n" );
	fprintf(stderr, "Octree[1]\n");
	printOctree_NB(trees[1]);
	fprintf(stderr, "---------------\n" );
	fprintf(stderr, "---------------\n" );

	//reduce
	_reduceOctrees_NB(threadPool, trees, NBODY_NPROCS);
	computeMassVals_NB(trees[0]->root);


	return retVals[0];
}


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
	long numN[NBODY_NPROCS])
{

	ExecutorThreadPool& threadPool = ExecutorThreadPool::getThreadPool();

	//map
	for (int i = 0; i < NBODY_NPROCS; ++i) {
		std::function<void()> f = std::bind(computeForcesOctreeBH_NB,
			numN[i],
			m + startN[i],
			r + 3*startN[i],
			work + startN[i],
			a + 3*startN[i],
			tree,
			list1[i],
			list2[i],
			thetaMAC
		);
		threadPool.addTaskAtIdx(f, i);
	}

	threadPool.waitForAllThreads(); //sync before returning;

}


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
	long numN[NBODY_NPROCS])
{

	ExecutorThreadPool& threadPool = ExecutorThreadPool::getThreadPool();

	//map
	for (int i = 0; i < NBODY_NPROCS; ++i) {
		std::function<void()> f = std::bind(performNBodyHalfStepA,
			numN[i],
			dt,
			r + 3*startN[i],
			v + 3*startN[i],
			a + 3*startN[i],
			m + startN[i]
		);
		threadPool.addTaskAtIdx(f, i);
	}

	threadPool.waitForAllThreads(); //sync before returning;

}


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
	long numN[NBODY_NPROCS])
{

	ExecutorThreadPool& threadPool = ExecutorThreadPool::getThreadPool();

	//map
	for (int i = 0; i < NBODY_NPROCS; ++i) {
		std::function<void()> f = std::bind(performNBodyHalfStepB,
			numN[i],
			dt,
			r + 3*startN[i],
			v + 3*startN[i],
			a + 3*startN[i],
			m + startN[i]
		);
		threadPool.addTaskAtIdx(f, i);
	}

	threadPool.waitForAllThreads(); //sync before returning;

}
