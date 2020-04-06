

#include "NBodyForces.h"
#include <math.h>

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
	double* __restrict__ a)
{
    //         a[k] = -r[k] / (r2*sqrtr2);
	double ai[3] = {0.0, 0.0, 0.0};
	double rx = r[3*i + 0];
	double ry = r[3*i + 1];
	double rz = r[3*i + 2];
	double dx, dy, dz, D;
	long j;
	for (j = 0; j < i; ++j) {
		//really dx is other way around, be this way we can avoid -1.0* later.
		dx = r[3*j] - rx;
		dy = r[3*j + 1] - ry;
		dz = r[3*j + 2] - rz;
		D = dx*dx + dy*dy + dz*dz;
		D += _SOFTENING*_SOFTENING;
		D = 1.0 / (D*sqrt(D));
		ai[0] += m[j]*dx*D;
		ai[1] += m[j]*dy*D;
		ai[2] += m[j]*dz*D;
	}
	for (j = i+1; j < n; ++j) {
		dx = r[3*j] - rx;
		dy = r[3*j + 1] - ry;
		dz = r[3*j + 2] - rz;
		D = dx*dx + dy*dy + dz*dz;
		D += _SOFTENING*_SOFTENING;
		D = 1.0 / (D*sqrt(D));
		ai[0] += m[j]*dx*D;
		ai[1] += m[j]*dy*D;
		ai[2] += m[j]*dz*D;
	}

	a[3*i + 0] = ai[0];
	a[3*i + 1] = ai[1];
	a[3*i + 2] = ai[2];
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
void computeForceMonoInteractList_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a,
	const NBOctreeNode_t** intList,
	long listSize)
{
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;

	const NBOctreeNode_t* node;
	double dx, dy, dz, D, D2;
	for(int i = 0; i < listSize; ++i) {
		node = intList[i];

		//convert r to c.o.m. referece frame of the node.
		dx = node->com[0] - r[0];
		dy = node->com[1] - r[1];
		dz = node->com[2] - r[2];
		// m*r / |r|^3 : normal monopole part/direct interaction
		D2 = dx*dx + dy*dy + dz*dz;
		D2 += _SOFTENING*_SOFTENING;
		D = 1.0 / sqrt(D2); //hopefully compiler catches the one over sqrt
		D2 = 1.0 / D2;
		D = D*D2; // 1/D3

		//m*r / |r|^3
		a[0] += node->mass*dx*D;
		a[1] += node->mass*dy*D;
		a[2] += node->mass*dz*D;

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
void computeForceInteractList_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a,
	const NBOctreeNode_t** intList,
	long listSize)
{
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;

	const NBOctreeNode_t* node;
	double dx, dy, dz, D, D2;
	double qx, qy, qz;
	for(int i = 0; i < listSize; ++i) {
		node = intList[i];

		//convert r to c.o.m. referece frame of the node.
		dx = node->com[0] - r[0];
		dy = node->com[1] - r[1];
		dz = node->com[2] - r[2];
		// m*r / |r|^3 : normal monopole part/direct interaction
		D2 = dx*dx + dy*dy + dz*dz;
		D2 += _SOFTENING*_SOFTENING;
		D = 1.0 / sqrt(D2); //hopefully compiler catches the one over sqrt
		D2 = 1.0 / D2;
		D = D*D2; // 1/D3

		//m*r / |r|^3
		a[0] += node->mass*dx*D;
		a[1] += node->mass*dy*D;
		a[2] += node->mass*dz*D;

		if (node->N > 1) {
			//just did monopole so now quadrupole approximate

			//Q.r / |r|^5; recall quadMom is only the upper triangle of symmetric tensor
			qx = node->quadMom[0]*dx + node->quadMom[1]*dy + node->quadMom[2]*dz;
			qy = node->quadMom[1]*dx + node->quadMom[3]*dy + node->quadMom[4]*dz;
			qz = node->quadMom[2]*dx + node->quadMom[4]*dy + node->quadMom[5]*dz;


			D *= D2; // 1/D5 now
			a[0] -= qx*D;
			a[1] -= qy*D;
			a[2] -= qz*D;

			//5*r.Q.r*r / 2*|r|^7
			qx = dx*qx + dy*qy + dz*qz;
			qx *= 2.5;
			D *= D2; // 1/D7 now

			a[0] += qx*dx*D;
			a[1] += qx*dy*D;
			a[2] += qx*dz*D;
		}
	}
}
