


#ifndef _NBODY_HASHOCTREE_
#define _NBODY_HASHOCTREE_

#include "NBodyConfig.h"
#include "NBodyHelpers.h"
#include "NBodyKeys.h"
#include "NBodyOctree.h"

#include <stdint.h>

#if NBODY_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


#define getOctParentKey_NB(childKey) ((childKey) >> 3)

#define getOctChildKey_NB(parentKey, OCT) ((parentKey) << 3 | (OCT))

#define getOctantFromKey_NB(key) ((key) & 7);

#define ROOT_KEY 1

static inline int isKeyContainedIn_NB(spatialKey_t key, spatialKey_t parent) {
	if (parent == 1) {
		return 1;
	}
	if (key < parent) {
		return 0;
	}
	int parLev = 0;
	int keyLev = 0;
	spatialKey_t p = parent;
	while ( p >>= 3 ) { parLev += 3; }
	p = key;
	while ( p >>= 3 ) { keyLev += 3; }
	key >>= (keyLev - parLev);
	return (key == parent);
}


/**
 * A node of the hashed octree for graviational N-body.
 */
typedef struct NBodyHOTNode {
	//extents of the node is, component-wise, [(center - size), (center + size)]
	double center[3];
	double size;

	//If N = 1 mass and com are in fact the single particle's mass and pos.
	double com[3];
	double mass;
	long N; //number of bodies at or below this node.

	//Upper triangle of quadrupole moment tensor Q ordered as
	//Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
	double quadMom[6];

	spatialKey_t key; //the spatial key for this node/body
	uint8_t childByte; //a bitfield encoding which children actually exist
	uint8_t childrenLocal;
	uint8_t requested;
	int owningRank; //the rank of the processor which owns this node

	struct NBodyHOTNode* next; // a pointer to the next node in case of collisions

} NBodyHOTNode_t;


/**
 * An Hashed Octree (liner octree) of NBodyHOTNodes.
 * Implemented using a hash table and per-node keys of the tree.
 * See "Linear Quadtree" in "An Effective Way to Represent Quadtrees" (Gargantini, 1982)
 * and "A Parallel Hashed Oct-Tree N-Body Algorithm" (Warren and Salmon, 1993).
 */
typedef struct NBodyHashedOctree {

	NBodyHOTNode_t** data;
	int h;
	int hash;

} NBodyHOT_t;


/**
 * Multipole acceptability criteria pointer.
 * Based on an octree node, a test point's position,
 * and an MAC parameter, return true iff criteria passes.
 * That is so say, if the multipole approximation provided by
 * the node is acceptable for the test point.
 */
typedef int (*HOTMAC_ptr)(const NBodyHOTNode_t*, const double*, double);


#if NBODY_MPI
/**
 * Struct to hold information required for asynchronous
 * data transfer during distributed tree traversal.
 */
typedef struct deferData {

	NBodyHOTNode_t* parentNode; //node being deferred and requested
	uint8_t* recvBuffer; //pointer posted async receive buffer
	MPI_Request reqs[2]; //MPI_request handles for the key request and incoming data

} NBodyDeferData_t;
#endif


/**
 * Create a hashed octree (HOT) with allocation 2^len and thus
 * a depth of len-3.
 *
 * @param len: 2^len is the size of the hash table in the HOT.
 * @return the newly created HOT.
 */
static NBodyHOT_t createHashedOctree_NB(int len) {
	NBodyHOT_t ret;
	ret.h = len;
	ret.hash = (1 << len) - 1;
	ret.data = (NBodyHOTNode_t**) calloc((1<<len), sizeof(NBodyHOTNode_t*));
	return ret;
}

/**
 * Free all memory allocated to a hashed octree.
 * This includes the hash table itself and all nodes contained in the table.
 *
 * @param hot: the hashed octree
 */
static inline void freeHashedOctree_NB(NBodyHOT_t hot) {
	if (hot.data != NULL) {
		int len = 1 << hot.h;
		NBodyHOTNode_t *node, *nextNode;
		for (int i = 0; i < len; ++i) {
			node = hot.data[i];
			//traverse chain to free it
			while (node != NULL) {
				nextNode = node->next;
				free(node);
				node = nextNode;
			}
		}
		free(hot.data);
	}
	hot.data = NULL; //to be safe
}



/**
 * Empty the contents of the given hashed octree,
 * freeing all nodes but retaining the hash table itself.
 *
 * @param hot: the hashed octree whose contents should be cleared
 */
static inline void clearHashedOctree_NB(NBodyHOT_t hot) {
	if (hot.data != NULL) {
		int len = 1 << hot.h;
		NBodyHOTNode_t *node, *nextNode;
		for (int i = 0; i < len; ++i) {
			node = hot.data[i];
			//traverse chain to free it
			while (node != NULL) {
				nextNode = node->next;
				free(node);
				node = nextNode;
			}
			hot.data[i] = NULL;
		}
	}
}


/**
 * Find a particular node/body in the hashed octree by key.
 *
 * @param hot: the hashed octree
 * @param key: the key of the node or body
 * @return the node, if found, else NULL.
 */
static NBodyHOTNode_t* getHOTNode_NB(NBodyHOT_t hot, spatialKey_t key) {
	if (hot.data == NULL) {
		return NULL;
	}

	int k = key & hot.hash;
	NBodyHOTNode_t* node = hot.data[k];
	while (node != NULL && node->key != key) {
		node = node->next;
	}
	return node;
}

/**
 * Insert a hashed octree node into the HOT hash table.
 *
 * @param hot: the hashed octree
 * @param node: the node to insert
 */
static void insertHOTNode_NB(NBodyHOT_t hot, NBodyHOTNode_t* node) {
	if (hot.data == NULL) {
		return;
	}

	int k = node->key & hot.hash;
	node->next = hot.data[k];
	hot.data[k] = node;
}

/**
 * Remove a hased octree node from the HOT hash table.
 *
 * @param hot: the hashed octree
 * @param key: the key of the node to remove
 * @return the removed hashed octree node
 */
static NBodyHOTNode_t* removeHOTNode_NB(NBodyHOT_t hot, spatialKey_t key) {
	if (hot.data == NULL) {
		return NULL;
	}

	int k = key & hot.hash;
	NBodyHOTNode_t* node = hot.data[k];
	NBodyHOTNode_t* prevNode = NULL;
	while (node->key != key) {
		prevNode = node;
		node = node->next;
	}
	if (prevNode != NULL) {
		prevNode->next = node->next;
	} else {
		hot.data[k] = node->next;
	}

	node->next = NULL;
	return node;
}

/**
 * Recursively free this node and all of its children.
 *
 * @param node: the root node of the (sub-)tree to free.
 */
static inline void freeHOTNode_NB(NBodyHOTNode_t* node) {
	free(node);
}

/**
 * Create an empty HOT node with default values.
 *
 * @return the newly constructed node.
 */
NBodyHOTNode_t* createHOTNode_NB();


/**
 * Create a child node of octree node with the specified parent.
 *
 * @param parent: the parent node in the hashed octree node.
 * @param oct: the child octant of parent the new node should inhabit.
 * @param halfSize: half the size of the parent node; the size of the child node.
 * @param r: a pointer to the triple giving the child node's body's position.
 * @param m: the mass of the body contained in this child node.
 * @return the newly create hashed octree node.
 */
NBodyHOTNode_t* createHOTChildNode_NB(NBodyHOTNode_t* parent, NBOctant_t oct,
		double halfSize, const double* r, double m);


/**
 * Given an octree node with exactly one point in it,
 * create a child node and push that point down to that child.
 * This is done in preparation of including another point.
 *
 * @param pNode: the node with one point to push to a child octant
 */
static inline void pushHOTNodePointToChild_NB(NBodyHOT_t hot, NBodyHOTNode_t* pNode) {
	if ( pNode == NULL || pNode->N != 1) {
		return;
	}

	//Since N == 1, com is the position of that one body.
	NBOctant_t oct = determineOctant_NB(pNode->center, pNode->com);
	NBodyHOTNode_t* child = createHOTChildNode_NB(pNode, oct, pNode->size*0.5, pNode->com, pNode->mass);
	child->owningRank = pNode->owningRank;
	setOctChild_NB(pNode->childByte, oct);
	insertHOTNode_NB(hot, child);
}


/**
 * Internal method for updating node's COM and multipole moments.
 */
void _computeMassValsHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* node);


/**
 * Given a filled in hashed octree,
 * compute mass, center of mass, and multipole moments
 * for each node.
 *
 * @param hot: the hashed octree
 */
static inline void computeMassValsHOT_NB(NBodyHOT_t hot) {
	_computeMassValsHOT_NB(hot, getHOTNode_NB(hot, ROOT_KEY));
}


/**
 * Given an hashed octree, with at least a root node already inserted,
 * insert the body at position r into the octree, splitting nodes as needed.
 *
 * @param hot: the hashed octree
 * @param node: root of subtree into which to insert, or NULL for root
 * @param r: the position of the body to insert.
 * @param m: the mass of the body to insert.
 */
void insertBodyHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, const double* r, double m);


/**
 * Build an octree "in place", re-using as much of the pre-allocated space as possible,
 * thus destroying the existing octree in the process.
 * Given a list of N positions and masses, bounded by (-domainSize, domainSize)
 * in each dimension, empty the octree pointed to tree_ptr and then
 * fill it with the N bodies.
 * If tree_ptr points to NULL, a new octree is created.
 * If tree_ptr points to a tree whose size is not domainSize, the existing tree
 * is freed and a new tree is created and returned in tree_ptr.
 *
 * @param N: the number of bodies.
 * @param r: an array of 3*N values representing the bodies' positions.
 * @param m: an array of N values representing the bodies' masses.
 * @param domainSize: the size of the entire octree.
 * @param tree_ptr: a pointer to the existing tree to update.
 * @return the new hashed octree
 */
NBodyHOT_t buildHashedOctreeInPlace_NB(long N, const double* r, const double* m, double domainSize, NBodyHOT_t hot);


/**
 * Given a list of N positions and masses, bounded by (-domainSize, domainSize)
 * in each dimension, create an octree containing all bodies.
 *
 * @param N: the number of bodies.
 * @param r: an array of 3*N values representing the bodies' positions.
 * @param m: an array of N values representing the bodies' masses.
 * @param domainSize: the size of the entire octree.
 * @param h: the allocation of the hashed octree as 2^h.
 * @return the newly created octree.
 */
static inline NBodyHOT_t buildHashedOctree_NB(long N, const double* r, const double* m, double domainSize, int h) {
	NBodyHOT_t hot = createHashedOctree_NB(h);
	return buildHashedOctreeInPlace_NB(N, r, m, domainSize, hot);
}


/**
 * Walk the octree start at root, determining the interactions
 * that should be including for the test point as position r.
 * This method works iteratively, without recursive calls.
 *
 * @param root: the roto of the tree to traverse.
 * @param r: an array of size 3 representing the position of the test point.
 * @param mac: a pointer to the multipole acceptance function.
 * @param theta: a parameter to the MAC.
 * @param walkList: an array of size at least root->N used as working space.
 * @param[out] interactList: an array of size at least root->N in which to
 *                           store the interacting cells.
 * @return the length of the interactList.
 */
long traverseHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	const NBodyHOTNode_t** walkList,
	const NBodyHOTNode_t** interactList);


/**
 * Merge two octrees in place, moving and inserting from srcTree to dstTree.
 * srcTree is emptied and free'd by there operation and should not be used afterward.
 * dstTree and srcTree must have same hash table size.
 *
 * @param dstTree, a pointer to the destination tree for the merge result
 * @param srcTree, the tree to merge cells from.
 */
void mergeHashedOctrees_NB(NBodyHOT_t* dstTree, NBodyHOT_t srcTree);


#if NBODY_MPI

/**
 * At the borders of each rank's domain (i.e. bodies with index 0 and localN-1),
 * exchange those bodies with the respective neighbour (assuming a linear topology)
 * so that each rank's octree is consistent in the placing
 * of its leaf nodes, as if global knowledge of the bodies was known.
 *
 * This process involves exchanging bodies on the border(s) of the rank's domain
 * and then inserting each neighbour's body into the local hashed octree.
 *
 * @param hot: the local processor's octree.
 * @param localN: the number of bodies local to this process.
 * @param r: the position data for all local bodies, of size 3*localN.
 * @param m: the mass data for all local bodies, of size localN.
 */
void exchangeBorderBodies_NB(NBodyHOT_t hot, long localN, double* r, double* m);


/**
 * Given a hashed octree which contains the border bodies of adjacent
 * processors, determine which other nodes are the "branch nodes",
 * the set of cells which represent the entire local hashed octree
 * at the coarsest level possible.
 * The branch nodes thus exclude the parts of the tree which contain
 * border bodies.
 *
 * @see exchangeBorderBodies_NB
 *
 * @param hot: the hashed octree
 * @return a linked list of the branch nodes
 *
 */
NBodyHOTNode_t* getBranchNodes_NB(NBodyHOT_t hot);


/**
 * Broadcast the branch nodes of each rank's hashed octree
 * to all other ranks, merging those branch nodes into the local hashed octree.
 *
 * @param hot: the local hashed octree
 * @param branchData: a pointer to a possibly preallocated array of bytes for serialization
 * @param branchAlloc: a pointer to the allocation size of the data array
 */
void exchangeBranchNodes_NB(NBodyHOT_t hot,	uint8_t** branchData, size_t* branchAlloc);


/**
 * In a distributed hashed octree, request child data
 * of a (child of a) branch node from another processor.
 *
 * @param hot: the local hashed octree
 * @param treeData: a pointer to a byte array for temporary serialization storage
 * @param treeAlloc: a pointer to the current allocation size of the byte array
 */
void requestAndInsertHOTChildren_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, uint8_t** treeData, size_t* treeAlloc);


/**
 * In a distributed hashed octree, check for, and respond to, requests
 * for child data of a (child of a) branch node.
 *
 * @param hot: the local hashed octree
 * @param treeData: a pointer to a byte array for temporary serialization storage
 * @param treeAlloc: a pointer to the current allocation size of the byte array
 */
void fulfillHOTChildrenRequests_NB(NBodyHOT_t hot, uint8_t** treeData, size_t* treeAlloc);


/**
 * Walk the octree starting at the root, determining the interactions
 * that should be including for the test point as position r.
 * The hashed octree is not necessarily fully present on the local processor,
 * requests for non-local data may be required.
 *
 * @param root: the root of the tree to traverse.
 * @param r: an array of size 3 representing the position of the test point.
 * @param mac: a pointer to the multipole acceptance function.
 * @param theta: a parameter to the MAC.
 * @param walkList: an array of size at least root->N used as working space.
 * @param[out] interactList: an array of size at least root->N in which to
 *                           store the interacting cells.
 * @return the length of the interactList.
 */
long traverseDistributedHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	NBodyHOTNode_t** walkList,
	NBodyHOTNode_t** interactList,
	uint8_t** treeData,
	size_t* treeAlloc);


long traverseDistribAsyncHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	NBodyHOTNode_t** walkList,
	// NBodyHOTNode_t** deferList,
	NBodyHOTNode_t** interactList,
	uint8_t** treeData,
	size_t* treeAlloc);
#endif


/**
 * Recrusively print to stderr the sub-octree whose root is at node.
 *
 * @param hot: the hashed octree encompassing all nodes
 * @param node: the root of the sub-octree to print.
 * @param level: the current level in the octree.
 */
static inline void printHOTNode_NB(NBodyHOT_t hot, const NBodyHOTNode_t* node, int level) {
	if (node == NULL) {
		return;
	}
	for(int i = 0; i < level; ++i) {
		fprintf(stderr, "--");
	}

	if (node->N == 1) {
		spatialKey_t bodyKey;
		computeSpatialKey_NB(0, node->com, &bodyKey, getHOTNode_NB(hot, ROOT_KEY)->size);
		fprintf(stderr, "(%.5f, %.5f, %.5f), %llo, %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", node->center[0], node->center[1], node->center[2], node->key, bodyKey, node->N, node->com[0], node->com[1], node->com[2], node->mass, node->owningRank);

	} else {
		fprintf(stderr, "(%.5f, %.5f, %.5f), %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", node->center[0], node->center[1], node->center[2], node->key, node->N, node->com[0], node->com[1], node->com[2], node->mass, node->owningRank);

	}


	for (int i = 0; i < 8; ++i) {
		if (hasOctChild_NB(node->childByte, i)) {
			printHOTNode_NB(hot, getHOTNode_NB(hot, getOctChildKey_NB(node->key, i)), level+1);
		}
	}
}


/**
 * Print a hashed octree to stderr using an in-order traversal.
 * @param hot: the hased octree to print.
 */
static inline void printHashedOctree_NB(NBodyHOT_t hot) {
	printHOTNode_NB(hot, getHOTNode_NB(hot, ROOT_KEY), 0);
}


/**
 * Determine the number of nodes in an octree.
 * @param tree: the octree.
 * @return the number of nodes in the tree.
 */
static long countHOTNodes_NB(NBodyHOT_t hot) {
	long nnodes = 0;
	if (hot.data != NULL) {
		int len = 1 << hot.h;
		NBodyHOTNode_t* node;
		for (int i = 0; i < len; ++i) {
			node = hot.data[i];
			while (node != NULL) {
				++nnodes;
				node = node->next;
			}
		}
	}
	return nnodes;
}



#ifdef __cplusplus
}
#endif


#endif
