
#ifndef _NBODY_FORCES_
#define _NBODY_FORCES_

#include <time.h>
#include <stdio.h>

#include "NBodyOctree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * A softening parameter for the gravitational potential calculation
 * to avoid huge fluctuations caused by close encounters.
 */
static const double _SOFTENING = 0.025;


/**
 * Compute the acceleration of body i from the mutual graviational
 * potential of the other n-1 bodies.
 *
 * @param n, the number of bodies
 * @param i, the index of the body of which to compute acceleration.
 * @param m, an array of n doubles holding the masses of the bodies
 * @param r, an array of 3*n doubles holding the positions of the bodies
 * @param a, an array of 3*n doubles to hold the resulting acceleration of the bodies.
 */
void computeForce(long n, long i,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a);


/**
 * Using gravitational potential, compute and update the accelerations
 * of the n bodies in place.
 *
 * @param n, the number of bodies
 * @param m, an array of n doubles holding the masses of the bodies
 * @param r, an array of 3*n doubles holding the positions of the bodies
 * @param[out] a, an array of 3*n doubles to hold the resulting acceleration of the bodies.
 */
static inline void computeForces(long n,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a) {
	for (long i = 0; i < n; ++i) {
		computeForce(n, i, m, r, a);
	}
}


/**
 * Determine is the multipole approximation is acceptable for the given node
 * and a test particle at position r based on the Barnes-Hut criteria
 * of opening angle theta.
 *
 * @param node, the octree node to compare against.
 * @param r, the position of the target particle to test.
 * @param theta, a parameter controling the MAC tolerance.
 * @return 1 iff the multipole approxiation is acceptable.
 */
static inline int barnesHutMAC_NB(const NBOctreeNode_t* node, const double* r, double theta) {
	// node's diamater / distance-to-point < theta
	// fprintf(stderr, "l = %.5f\n", 2.0*node->size);
	// fprintf(Stderr, "d = %.5f\n", )
	float dx = node->com[0]-r[0];
	float dy = node->com[1]-r[1];
	float dz = node->com[2]-r[2];
	dx = dx*dx + dy*dy + dz*dz;
	return (4.0*node->size*node->size < theta*theta*dx);
}


/**
 * Compute the acceleration for a particle at position r with mass m
 * based on the list of Octree nodes in interaction list.
 * This nodes in the interaction list may be leaf nodes or internal nodes,
 * as determined by a MAC.
 * If a node is internal, the quadrupole moments may be included in the potential
 * calculation.
 *
 * @param m, a pointer to the mass of the body whose acceleration is being computes.
 * @param r, the 3-dim position of the body.
 * @param[out] a, an array of size at least 3 to store the computed acceleration.
 * @param intList, the list of interacting octree nodes.
 * @param listSize, the size of the interaction list.
 */
void computeForceInteractList_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a,
	const NBOctreeNode_t** intList,
	long listSize);


/**
 * Compute the acceleration for n bodies at positions r with mass m
 * making use of the octree and a Barnes-Hut MAC for multipole approximation.
 * If a node is internal, the quadrupole moments are included in the potential
 * calculation.
 *
 * @param n, the number of bodies
 * @param m, an array of n doubles holding the masses of the bodies
 * @param r, an array of 3*n doubles holding the positions of the bodies
 * @param[out] a, an array of 3*n doubles to hold the resulting acceleration of the bodies.
 * @param tree, an octree holding the n bodies.
 * @param list1, an array used as working space of size at least N.
 * @param list2, an array used as working space of size at least N.
 * @param thetaMac, a MAC parameter to control the use of direct calculation or multipole approximation.
 */
static inline void computeForcesOctreeBH_NB(long n,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* work,
	double* __restrict__ a,
	const NBOctree_t* tree,
	const NBOctreeNode_t** list1,
	const NBOctreeNode_t** list2,
	double thetaMAC)
{
	long listLen;
	for (long i = 0; i < n; ++i) {
		listLen = traverseTreeInteractionList_NB(tree, r + 3*i, &(barnesHutMAC_NB), thetaMAC, list1, list2);
		work[i] = 0.0;
		for (long j = 0; j < listLen; ++j) {
			//since quadrupole moments double flops for interaction, weight those double.
			work[i] += ((list2[j]->N >> 1) & 0x1) + 1; // (N > 1) + 1
		}
		// listLen = traverseTreeInteractionList_NB(tree, r + 3*i, &(barnesHutMAC_NB), thetaMAC, list1, list2);
		computeForceInteractList_NB(m + i, r + 3*i, a + 3*i, list2, listLen);
	}
}


/**
 * Compute the acceleration for a particle at position r with mass m
 * based on the list of Octree nodes in interaction list.
 * This nodes in the interaction list may be leaf nodes or internal nodes,
 * as determined by a MAC.
 * If a node is internal, only the monopole moments are included in the potential
 * calculation.
 *
 * @param m, a pointer to the mass of the body whose acceleration is being computes.
 * @param r, the 3-dim position of the body.
 * @param[out] a, an array of size at least 3 to store the computed acceleration.
 * @param intList, the list of interacting octree nodes.
 * @param listSize, the size of the interaction list.
 */
void computeForceMonoInteractList_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a,
	const NBOctreeNode_t** intList,
	long listSize);


/**
 * Compute the acceleration for n bodies at positions r with mass m
 * making use of the octree and a Barnes-Hut MAC for multipole approximation.
 * If a node is internal, only the monopole moments are included in the potential
 * calculation.
 *
 * @param n, the number of bodies
 * @param m, an array of n doubles holding the masses of the bodies
 * @param r, an array of 3*n doubles holding the positions of the bodies
 * @param[out] a, an array of 3*n doubles to hold the resulting acceleration of the bodies.
 * @param tree, an octree holding the n bodies.
 * @param list1, an array used as working space of size at least N.
 * @param list2, an array used as working space of size at least N.
 * @param thetaMac, a MAC parameter to control the use of direct calculation or multipole approximation.
 */
static inline void computeForcesMonoOctreeBH_NB(long n,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* work,
	double* __restrict__ a,
	const NBOctree_t* tree,
	const NBOctreeNode_t** list1,
	const NBOctreeNode_t** list2,
	double thetaMAC)
{
	long listLen;
	for (long i = 0; i < n; ++i) {
		listLen = traverseTreeInteractionList_NB(tree, r + 3*i, &(barnesHutMAC_NB), thetaMAC, list1, list2);
		work[i] = (double) listLen;
		computeForceMonoInteractList_NB(m + i, r + 3*i, a + 3*i, list2, listLen);
	}
}

#ifdef __cplusplus
}
#endif

#endif
