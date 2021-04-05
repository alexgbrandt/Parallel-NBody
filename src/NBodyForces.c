

#include "NBodyForces.h"
#include <math.h>

#if NBODY_MPI
#include <mpi.h>
#include <string.h>
#endif

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
	double* a)
{
    //         a[k] = -r[k] / (r2*sqrtr2);
	double* __restrict__ ai = a + 3*i;
	ai[0] = 0.0; ai[1] = 0.0; ai[2] = 0.0;
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


/**
 * Compute the acceleration for a particle at position r with mass m
 * based on the list of hased octree nodes in interaction list.
 * This nodes in the interaction list may be leaf nodes or internal nodes,
 * as determined by a MAC.
 * If a node is internal, only the monopole moments are included in the potential
 * calculation.
 *
 * @param m, a pointer to the mass of the body whose acceleration is being computes.
 * @param r, the 3-dim position of the body.
 * @param[out] a, an array of size at least 3 to store the computed acceleration.
 * @param intList, the list of interacting hased octree nodes.
 * @param listSize, the size of the interaction list.
 */
void computeForceHOTInteractList_NB(
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* __restrict__ a,
	const NBodyHOTNode_t** intList,
	long listSize)
{
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;

	const NBodyHOTNode_t* node;
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




#if NBODY_MPI
/**
 * Update the accelerations of bodies from the mutual graviational
 * potential of the other n-1 bodies.
 *
 * @param localN: the number of bodies on each processor
 * @param m: an array of n doubles holding the masses of the local bodies
 * @param r: an array of 3*n doubles holding the positions of the local bodies
 * @param a: an array of 3*n doubles to hold the resulting acceleration of the local bodies.
 * @param m_rec: array of n doubles each to use as a message buffer
 * @param r_rec: array of 3*n doubles each to use as a message buffer
 */
void computeForcesNaive_NBMPI(long localN,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* a,
	double* __restrict__ m_rec,
	double* __restrict__ r_rec)
{

	//first update accels by local bodies
	// a[k] = -r[k] / (r2*sqrtr2);

	double* __restrict__ ai;
	double rx, ry, rz;
	long i, j;
	double dx, dy, dz, D;

	for (i = 0; i < localN; ++i) {
		rx = r[3*i + 0];
		ry = r[3*i + 1];
		rz = r[3*i + 2];
		ai = a + 3*i;
		ai[0] = 0.0; ai[1] = 0.0; ai[2] = 0.0;

		for (j = 0; j < i; ++j) {
			//really dx is other way around, but this way we can avoid -1.0* later.
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
		for (j = i+1; j < localN; ++j) {
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
	}

	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int src_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;
	int dst_rank = (world_rank + 1) % world_size;
	MPI_Status stat;

	//fill send buffers initially
	memcpy(m_rec, m, sizeof(double)*localN);
	memcpy(r_rec, r, sizeof(double)*localN*3);

	for (int rounds = 1; rounds < world_size; ++rounds) {
		// fprintf(stderr, "Current rank: %d, src rank: %d, dest rank: %d\n", world_rank, dst_rank, src_rank);
		MPI_Sendrecv_replace(m_rec, localN, MPI_DOUBLE,
							 dst_rank, 0,
							 src_rank, 0,
							 MPI_COMM_WORLD, &stat);
		MPI_Sendrecv_replace(r_rec, localN*3, MPI_DOUBLE,
							 dst_rank, 0,
							 src_rank, 0,
							 MPI_COMM_WORLD, &stat);

		//update each local body with the influence from the new bodies
		//TODO could use blocking here for locality
		for (i = 0; i < localN; ++i) {
			rx = r[3*i + 0];
			ry = r[3*i + 1];
			rz = r[3*i + 2];
			ai = a + 3*i;

			for (j = 0; j < localN; ++j) {
				//really dx is other way around, but this way we can avoid -1.0* later.
				dx = r_rec[3*j] - rx;
				dy = r_rec[3*j + 1] - ry;
				dz = r_rec[3*j + 2] - rz;
				D = dx*dx + dy*dy + dz*dz;
				D += _SOFTENING*_SOFTENING;
				D = 1.0 / (D*sqrt(D));
				ai[0] += m_rec[j]*dx*D;
				ai[1] += m_rec[j]*dy*D;
				ai[2] += m_rec[j]*dz*D;
			}
		}
	}
}

void computeForcesDistributedHOTBH_NB(long n,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* work,
	double* __restrict__ a,
	NBodyHOT_t hot,
	NBodyHOTNode_t** list1,
	NBodyHOTNode_t** list2,
	uint8_t** treeData,
	size_t* treeAlloc,
	double thetaMAC)
{
	long listLen;
	for (long i = 0; i < n; ++i) {
		listLen = traverseDistributedHOTInteractionList_NB(hot, r + 3*i, &(barnesHutHOTMAC_NB), thetaMAC,
					list1, list2, treeData, treeAlloc);

		if (work != NULL) {
			work[i] = 0.0;
			for (long j = 0; j < listLen; ++j) {
				//since quadrupole moments double flops for interaction, weight those double.
				work[i] += ((list2[j]->N >> 1) & 0x1) + 1; // (N > 1) + 1
			}
		}
		computeForceHOTInteractList_NB(m + i, r + 3*i, a + 3*i, (const NBodyHOTNode_t**) list2, listLen);

		//check for requested data
		fulfillHOTChildrenRequests_NB(hot, treeData, treeAlloc);
	}

	int allDone = 0;
	MPI_Request req;
	MPI_Ibarrier(MPI_COMM_WORLD, &req);
	while (!allDone) {
		fulfillHOTChildrenRequests_NB(hot, treeData, treeAlloc);
		MPI_Test(&req, &allDone, MPI_STATUS_IGNORE);
	}
}



void computeForcesAsyncDistributedHOTBH_NB(long n,
	const double* __restrict__ m,
	const double* __restrict__ r,
	double* work,
	double* __restrict__ a,
	NBodyHOT_t hot,
	NBodyHOTNode_t** list1,
	NBodyHOTNode_t** list2,
	uint8_t** treeData,
	size_t* treeAlloc,
	double thetaMAC)
{
	long listLen;
	for (long i = 0; i < n; ++i) {
		// listLen = traverseDistributedHOTInteractionList_NB(hot, r + 3*i, &(barnesHutHOTMAC_NB), thetaMAC,
		// 			list1, list2, treeData, treeAlloc);
		listLen = traverseDistribAsyncHOTInteractionList_NB(hot, r + 3*i, &(barnesHutHOTMAC_NB), thetaMAC,
					list1, list2, treeData, treeAlloc);

		if (work != NULL) {
			work[i] = 0.0;
			for (long j = 0; j < listLen; ++j) {
				//since quadrupole moments double flops for interaction, weight those double.
				work[i] += ((list2[j]->N >> 1) & 0x1) + 1; // (N > 1) + 1
			}
		}
		computeForceHOTInteractList_NB(m + i, r + 3*i, a + 3*i, (const NBodyHOTNode_t**) list2, listLen);

		//check for requested data
		asyncFulfillHOTChildrenRequests_NB(hot, *treeData);

	}

	int allDone = 0;
	MPI_Request req;
	MPI_Ibarrier(MPI_COMM_WORLD, &req);
	while (!allDone) {
		asyncFulfillHOTChildrenRequests_NB(hot, *treeData);
		MPI_Test(&req, &allDone, MPI_STATUS_IGNORE);
	}
}


#endif
