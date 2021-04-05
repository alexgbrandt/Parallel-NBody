


#ifndef _NBODY_OCTREE_
#define _NBODY_OCTREE_

#include "NBodyHelpers.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Define octree octants order to match morton ordering.
 */
#define XYZ_ORDER 1

/**
 * Enumerate the 8 octants, considering the parent at the origin
 * in an Open-GL oriented cooridate system:
 * (+X = right), (+Y = up), (+Z = out of the screen).
 * They can be used as symbolic indices into a node's children array.
 * Enumeration based on Morton ordering with X > Y > Z,
 * so that a triple (x,y,z) is interleaved ... | x_1y_1z_1 | x_0y_0z_0
 */
#if XYZ_ORDER
typedef enum NBodyOctant {
	NBOct_0 = 0x0, //---
	NBOct_1 = 0x1, //--+
	NBOct_2 = 0x2, //-+-
	NBOct_3 = 0x3, //-++
	NBOct_4 = 0x4, //+--
	NBOct_5 = 0x5, //+-+
	NBOct_6 = 0x6, //++-
	NBOct_7 = 0x7, //+++
} NBOctant_t;
#else
typedef enum NBodyOctant {
	NBOct_0 = 0x0, //+++
	NBOct_1 = 0x1, //-++
	NBOct_2 = 0x2, //--+
	NBOct_3 = 0x3, //+-+
	NBOct_4 = 0x4, //+--
	NBOct_5 = 0x5, //---
	NBOct_6 = 0x6, //-+-
	NBOct_7 = 0x7, //++-
} NBOctant_t;
#endif

/**
 * The direction in each dimension of each octant
 * w.r.t the center of the parent, indexed by
 * NBOctant_t.
 */
#if XYZ_ORDER
static const double OctantDir[] = {
	-1.0, -1.0, -1.0,
	-1.0, -1.0,  1.0,
	-1.0,  1.0, -1.0,
	-1.0,  1.0,  1.0,
	 1.0, -1.0, -1.0,
	 1.0, -1.0,  1.0,
	 1.0,  1.0, -1.0,
	 1.0,  1.0,  1.0
};
#else
static const double OctantDir[] = {
	 1.0,  1.0,  1.0,
	-1.0,  1.0,  1.0,
	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	 1.0, -1.0, -1.0,
	-1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0
};
#endif

#define hasOctChild_NB(CHILDBYTE, OCT) ((CHILDBYTE) & (1 << (OCT)))

#define setOctChild_NB(CHILDBYTE, OCT) ((CHILDBYTE) |= (1 << (OCT)))

#define unsetOctChild_NB(CHILDBYTE, OCT) ((CHILDBYTE) &= ~(1 << (OCT)))


/**
 * A node of an Octree for gravitational N-Body.
 *
 * Each node has a physical center and size.
 * Size is understood to be such that the octree node spans
 * (center[i] - size, center[i] + size), for each dimension i.
 *
 * Holds mass and center of mass over its eight children.
 */
typedef struct NBodyOctreeNode {
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

	uint8_t childByte;

	struct NBodyOctreeNode** children;

} NBOctreeNode_t;


/**
 * An Octree for gravitational N-Body.
 */
typedef struct NBodyOctree {

	NBOctreeNode_t* root;

} NBOctree_t;


/**
 * Multipole acceptability criteria pointer.
 * Based on an octree node, a test point's position,
 * and an MAC parameter, return true iff criteria passes.
 * That is so say, if the multipole approximation provided by
 * the node is acceptable for the test point.
 */
typedef int (*MAC_ptr)(const NBOctreeNode_t*, const double*, double);


/**
 * Create an empty octree node with default values.
 *
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyNode();


/**
 * Create an empty octree node whose size is the input size and whose center
 * is at the origin.
 *
 * @param size: the input size.
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyRootNode(double size);


/**
 * Create an empty octree node with a particular center and size.
 *
 * @param center: the node's center.
 * @param size: the node's size
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyEmptyNode(const double* center, double size);


/**
 * Given a parent octree node, create a particular child node given the octant
 * for that child and the particle's position and mass to be placed in that
 * child node.
 *
 * @param parent: the parent octree node
 * @param oct: the octant to place the new child
 * @param r: an array of size 3, the position of the body in the child.
 * @param m: the mass of the body in the child.
 * @return the newly created child.
 */
NBOctreeNode_t* createNBodyChildNode(
	NBOctreeNode_t* parent,
	NBOctant_t oct,
	const double* r,
	double m);


/**
 * Create an empty list to store octree child nodes.
 *
 * @return an array large enough to store octree node children.
 */
static inline NBOctreeNode_t** createNBodyChildList() {
	return (NBOctreeNode_t**) calloc(8, sizeof(NBOctreeNode_t*));
}


/**
 * Given a node's center, determine in which octant the position r belongs.
 *
 * @param center: the node's center.
 * @param r: the position to consider.
 * @return the octant which position r belongs.
 */
static inline NBOctant_t determineOctant_NB(double* center, const double* r) {
#if XYZ_ORDER
	int posx = (r[0] > center[0]);
	int posy = (r[1] > center[1]);
	int posz = (r[2] > center[2]);
	return  (NBOctant_t) ((posx << 2) | (posy << 1) | posz);
#else
	int negx = (r[0] - center[0]) < 0;
	int negy = (r[1] - center[1]) < 0;
	int negz = (r[2] - center[2]) < 0;
	int pos = (negx << 2) | (negy << 1) | negz;
	//pos as octant using greycode where - is a 1.
	if (pos == 0) {
		return NBOct_0;
	} else if(pos == 1) {
		return NBOct_7;
	} else if(pos == 2) {
		return NBOct_3;
	} else if(pos == 3) {
		return NBOct_4;
	} else if(pos == 4) {
		return NBOct_1;
	} else if(pos == 5) {
		return NBOct_6;
	} else if(pos == 6) {
		return NBOct_2;
	} else if(pos == 7) {
		return NBOct_5;
	}
	return NBOct_0;
#endif
}


/**
 * Given a parent octree node, and an octant, determine that octant's
 * center. Center is returned within octCenter
 *
 * @param parent: the parent octree node.
 * @param oct: the octant to consider.
 * @param[out] octCenter: an array of size 3 to hold the octant's center.
 */
static inline void determineOctantCenter_NB(
	NBOctreeNode_t* parent,
	NBOctant_t oct,
	double* octCenter)
{
	double halfSize = parent->size * 0.5;
	octCenter[0] = parent->center[0] + (halfSize*OctantDir[3*oct]);
	octCenter[1] = parent->center[1] + (halfSize*OctantDir[3*oct+1]);
	octCenter[2] = parent->center[2] + (halfSize*OctantDir[3*oct+2]);
}


/**
 * Given an octree node with exactly one point in it,
 * create a child node and push that point down to that child.
 * This is done in preparation of including another point.
 *
 * @param parent: the node with one point to push.
 */
static inline void pushNodePointToChild_NB(NBOctreeNode_t* parent) {
	if (parent == NULL || parent->N != 1) {
		return;
	}

	if (parent->children == NULL) {
		parent->children = createNBodyChildList();
	}
	NBOctant_t oct = determineOctant_NB(parent->center, parent->com);
	createNBodyChildNode(parent, oct, parent->com, parent->mass);
}


/**
 * Given a filled in octree whose root is node, recursively
 * compute mass, center of mass, and multipole moments
 * for each interior node and the root itself.
 *
 * @param node: the root of a (sub-)octree to compute values for.
 */
void computeMassVals_NB(NBOctreeNode_t* node);


/**
 * Given an (sub-)octree with root at node, insert
 * the body at position r into the octree, splitting nodes as needed.
 *
 * @param node: the root of the octree in which to insert the point.
 * @param r: the position of the body to insert.
 * @param m: the mass of the body to insert.
 */
void insertBodyOctree_NB(NBOctreeNode_t* node, const double* r, double m);

void insertNeighbourBody_NB(NBOctreeNode_t* node, const double* r, double m);

/**
 * Recursively free this node and all of its children.
 * @param node: the root node of the (sub-)tree to free.
 */
static inline void freeOctreeNode_NB(NBOctreeNode_t* node) {
	if (node != NULL) {
		if (node->children != NULL) {
			for(int i = 0; i < 8; ++i) {
				freeOctreeNode_NB(node->children[i]);
			}
			free(node->children);
		}
		free(node);
	}
}


/**
 * Free an Octree and all of its children recursively.
 * Actual freeing only occurs if reference count is low enough.
 *
 * @param tree: the octree to free.
 */
static inline void freeOctree_NB(NBOctree_t* tree) {
	if (tree != NULL) {
		freeOctreeNode_NB(tree->root);
		free(tree);
	}
}


/**
 * Given a (sub-)Octree, free all nodes which are empty.
 *
 * @param node: the root node of the (sub-)tree ot free.
 */
void pruneEmptyChildren_NB(NBOctreeNode_t* node);


/**
 * Empty the contents of the given octree node
 * and all of its children recursively.
 *
 * @param node: the node to empty
 */
static inline void emptyOctree_NB(NBOctreeNode_t* node) {
	if (node == NULL) {
		return;
	}
	node->N = 0;
	if (node->children != NULL) {
		for (int i = 0; i < 8; ++i) {
			emptyOctree_NB(node->children[i]);
		}
	}
}


/**
 * Given a list of N positions and masses, bounded by (-domainSize, domainSize)
 * in each dimension, create an octree containing all bodies.
 *
 * @param N: the number of bodies.
 * @param r: an array of 3*N values representing the bodies' positions.
 * @param m: an array of N values representing the bodies' masses.
 * @param domainSize: the size of the entire octree.
 * @return the newly created octree.
 */
NBOctree_t* buildOctree_NB(long N, const double* r, const double* m, double domainSize);


/**
 * Build an octree in place if possible.
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
 * @return 1 iff the build actually occurred in place.
 */
int buildOctreeInPlace_NB(long N, const double* r, const double* m, double domainSize, NBOctree_t** tree_ptr);


/**
 * Merge two sub-octrees with roots at dst and src.
 * Child nodes and moved and inserted from src to dst.
 * src is left with many or all nodes removed and moved to dst.
 * Both dst and src should have children already, in order to maintain
 * linkage from parent node to child node without adding reverse links.
 *
 * @param dst: the destination node for the merge result.
 * @param src: the source node to mrege nodes from.
 */
void mergeOctreeNodes_NB(NBOctreeNode_t* dst, NBOctreeNode_t* src);


/**
 * Merge two octrees in place, moving and inserting from srcTree to dstTree.
 * srcTree is left with many or all of its nodes removed and moved to dstTree.
 *
 * @param dstTree: the destination tree for the merge result
 * @param srcTree: the tree to merge cells from.
 */
void mergeOctreesInPlace_NB(NBOctree_t* dstTree, NBOctree_t* srcTree);


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
long traverseTreeInteractionList_NB(const NBOctree_t* root,
	const double* r,
	MAC_ptr mac,
	double theta,
	const NBOctreeNode_t** walkList,
	const NBOctreeNode_t** interactList);


/**
 * Walk the octree start at root, determining the interactions
 * that should be including for the test point as position r.
 * This method works recursively.
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
long traverseTreeInteractionListRec_NB(const NBOctree_t* tree,
	const double* r,
	MAC_ptr mac,
	double theta,
	const NBOctreeNode_t** walkList,
	const NBOctreeNode_t** interactList);


#if NBODY_MPI
/**
 * Using the map-reduce pattern, build an octree in parallel.
 * Rank 0 is the destination of the reduction and contains the entire tree.
 *
 * @param r: an array of 3*localN values for position.
 * @param m: an array of localN values for mass.
 * @param localN: the number of bodies to insert into the octree.
 * @param domainSize: the spatial size of the entire octree.
 * @param[in,out] tree: a pointer to the resulting octree. If it poitns to existing allocation, may be reused.
 * @param[in,out] treeData: message buffer for tree serialization, may change in size during this process
 * @param[in,out] treeSize: current allocation size of treeData, returns new treeData allocation
 *
 */
void mapReduceBuildOctrees_NBMPI(
	const double* __restrict__ r,
	const double* __restrict__ m,
	long localN,
	double domainSize,
	NBOctree_t** tree,
	uint8_t** treeData,
	size_t* treeSize);
#endif




/**
 * Recrusively print to stderr the sub-octree whose root is at node.
 *
 * @param node: the root of the sub-octree to print.
 * @param level: the current level in the octree.
 */
static inline void printOctreeNode_NB(const NBOctreeNode_t* node, int level) {
	if (node == NULL) {
		return;
	}
	for(int i = 0; i < level; ++i) {
		fprintf(stderr, "--");
	}
	fprintf(stderr, "(%.5f, %.5f, %.5f), N = %ld, COM=(%.5f, %.5f, %.5f) M=%.5f\n", node->center[0], node->center[1], node->center[2], node->N, node->com[0], node->com[1], node->com[2], node->mass);

	if (node->children != NULL) {
		for (int i = 0; i < 8; ++i) {
			printOctreeNode_NB(node->children[i], level+1);
		}
	}
}


/**
 * Print an octree to stderr.
 * @param tree: the tree to print.
 */
static inline void printOctree_NB(const NBOctree_t* tree) {
	printOctreeNode_NB(tree->root, 0);
}


/**
 * Count the number of levels in the octree rooted at node.
 * @param node: the root of the (sub-)octree.
 * @return the depth of the tree.
 */
static int countLevelsNode_NB(const NBOctreeNode_t* node) {
	if (node == NULL) {
		return 0;
	}
	if (node->children == NULL) {
		return 1;
	}

	int max = 0;
	for (int i = 0; i < 8; ++i) {
		max = MAX(max, countLevelsNode_NB(node->children[i]));
	}
	return max + 1;
}


/**
 * Determine the depth of the octree.
 * @param tree: the octree.
 * @return the depth of the tree.
 */
static int countLevels_NB(const NBOctree_t* tree) {
	if (tree == NULL) {
		return 0;
	}
	return countLevelsNode_NB(tree->root);
}


/**
 * Count the number of nodes in the octree rooted at node.
 * @param node: the root of the (sub-)octree.
 * @return the number of nodes in the tree.
 */
static long countNodesNode_NB(const NBOctreeNode_t* node) {
	if (node == NULL) {
		return 0;
	}
	if (node->children == NULL) {
		return 1;
	}

	int total = 1;
	for (int i = 0; i < 8; ++i) {
		total += countNodesNode_NB(node->children[i]);
	}
	return total;
}


/**
 * Determine the number of nodes in an octree.
 * @param tree: the octree.
 * @return the number of nodes in the tree.
 */
static long countNodes_NB(const NBOctree_t* tree) {
	if (tree == NULL) {
		return 0;
	}
	return countNodesNode_NB(tree->root);
}


#ifdef __cplusplus
}
#endif


#endif
