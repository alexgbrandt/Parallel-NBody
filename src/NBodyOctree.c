
#include "NBodyOctree.h"

#if NBODY_MPI
#include <mpi.h>
#include "NBodySerialize.h"
#endif

/**
 * Create an empty octree node with default values.
 *
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyNode() {
	NBOctreeNode_t* node = (NBOctreeNode_t*) malloc(sizeof(NBOctreeNode_t));
	node->center[0] = 0.0;
	node->center[1] = 0.0;
	node->center[2] = 0.0;
	node->size = 0.0;
	node->com[0] = 0.0;
	node->com[1] = 0.0;
	node->com[2] = 0.0;
	node->mass = 0.0;
	node->N = 0;

	node->quadMom[0] = 0.0;
	node->quadMom[1] = 0.0;
	node->quadMom[2] = 0.0;
	node->quadMom[3] = 0.0;
	node->quadMom[4] = 0.0;
	node->quadMom[5] = 0.0;

	node->childByte = 0;

	node->children = NULL;
	return node;
}


/**
 * Create an empty octree node whose size is the input size and whose center
 * is at the origin.
 *
 * @param size, the input size.
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyRootNode(double size) {
	NBOctreeNode_t* node = (NBOctreeNode_t*) malloc(sizeof(NBOctreeNode_t));
	node->center[0] = 0.0;
	node->center[1] = 0.0;
	node->center[2] = 0.0;
	node->size = size;
	node->com[0] = 0.0;
	node->com[1] = 0.0;
	node->com[2] = 0.0;
	node->mass = 0.0;
	node->N = 0;

	node->quadMom[0] = 0.0;
	node->quadMom[1] = 0.0;
	node->quadMom[2] = 0.0;
	node->quadMom[3] = 0.0;
	node->quadMom[4] = 0.0;
	node->quadMom[5] = 0.0;

	node->childByte = 0;

	node->children = NULL;
	return node;
}


/**
 * Create an empty octree node with a particular center and size.
 *
 * @param center, the node's center.
 * @param size, the node's size
 * @return the newly constructed node.
 */
NBOctreeNode_t* createNBodyEmptyNode(const double* center, double size) {
	NBOctreeNode_t* node = (NBOctreeNode_t*) malloc(sizeof(NBOctreeNode_t));
	node->center[0] = center[0];
	node->center[1] = center[1];
	node->center[2] = center[2];
	node->size = size;
	node->com[0] = 0.0;
	node->com[1] = 0.0;
	node->com[2] = 0.0;
	node->mass = 0.0;
	node->N = 0;

	node->quadMom[0] = 0.0;
	node->quadMom[1] = 0.0;
	node->quadMom[2] = 0.0;
	node->quadMom[3] = 0.0;
	node->quadMom[4] = 0.0;
	node->quadMom[5] = 0.0;

	node->childByte = 0;

	node->children = NULL;
	return node;
}


/**
 * Given a parent octree node, create a particular child node given the octant
 * for that child and the particle's position and mass to be placed in that
 * child node.
 *
 * @param parent, the parent octree node
 * @param oct, the octant to place the new child
 * @param r, an array of size 3, the position of the body in the child.
 * @param m, the mass of the body in the child.
 * @return the newly created child.
 */
NBOctreeNode_t* createNBodyChildNode(
		NBOctreeNode_t* parent,
		NBOctant_t oct,
		const double* r,
		double m)
{
	if (parent == NULL || parent->children[oct] != NULL) {
		return NULL;
	}
	double halfSize = parent->size * 0.5;
	parent->children[oct] = createNBodyNode();
	parent->children[oct]->center[0] = parent->center[0] + (halfSize*OctantDir[3*oct]);
	parent->children[oct]->center[1] = parent->center[1] + (halfSize*OctantDir[3*oct+1]);
	parent->children[oct]->center[2] = parent->center[2] + (halfSize*OctantDir[3*oct+2]);
	parent->children[oct]->size = halfSize;
	parent->children[oct]->com[0] = r[0];
	parent->children[oct]->com[1] = r[1];
	parent->children[oct]->com[2] = r[2];
	parent->children[oct]->mass = m;
	parent->children[oct]->N = 1;
	setOctChild_NB(parent->childByte, oct);
	return parent->children[oct];
}

/**
 * Given a filled in octree whose root is not, recursively
 * compute mass and conter of mass for each interior node and the root itself.
 *
 * @param node, the root of a (sub-)octree to compute values for.
 */
void computeMassVals_NB(NBOctreeNode_t* node) {
	if (node == NULL || node->children == NULL || node->N == 1) {
		return;
	}

	//we are now definitely an internal cell of the tree.

	//recursively update children first.
	int i;
	for (i = 0; i < 8; ++i) {
		computeMassVals_NB(node->children[i]);
	}

	//reset current node before accumulating
	node->mass = 0.0;
	node->com[0] = 0.0;
	node->com[1] = 0.0;
	node->com[2] = 0.0;

	//compute c.o.m. first.
	double mass;
	for (i = 0; i < 8; ++i) {
		if (node->children[i] != NULL) {
			mass = node->children[i]->mass;
			node->mass += mass;
			node->com[0] += mass*node->children[i]->com[0];
			node->com[1] += mass*node->children[i]->com[1];
			node->com[2] += mass*node->children[i]->com[2];

		}
	}
	node->com[0] /= node->mass;
	node->com[1] /= node->mass;
	node->com[2] /= node->mass;


	//compute quad moment contribution
	node->quadMom[0] = 0.0;
	node->quadMom[1] = 0.0;
	node->quadMom[2] = 0.0;
	node->quadMom[3] = 0.0;
	node->quadMom[4] = 0.0;
	node->quadMom[5] = 0.0;

	double pos[3];
	double xx, yy, zz, d2;
	for (i = 0; i < 8; ++i) {
		if (node->children[i] != NULL) {

			//first, get child's pos in c.o.m. reference frame
			pos[0] = node->children[i]->com[0] - node->com[0];
			pos[1] = node->children[i]->com[1] - node->com[1];
			pos[2] = node->children[i]->com[2] - node->com[2];

			mass = node->children[i]->mass;
			xx = pos[0]*pos[1]; // xy;
			node->quadMom[1] += 3.0*mass*xx;
			zz = pos[0]*pos[2]; // xz;
			node->quadMom[2] += 3.0*mass*zz;
			yy = pos[1]*pos[2]; // yz;
			node->quadMom[4] += 3.0*mass*yy;

			xx = pos[0]*pos[0];
			yy = pos[1]*pos[1];
			zz = pos[2]*pos[2];
			d2 = xx + yy + zz;
			d2 *= mass;
			mass *= 3.0;
			node->quadMom[0] += mass*xx - d2;
			node->quadMom[3] += mass*yy - d2;
			node->quadMom[5] += mass*zz - d2;

			//directly add child's moments also
			if (node->children[i]->N > 1) {
				node->quadMom[0] += node->children[i]->quadMom[0];
				node->quadMom[1] += node->children[i]->quadMom[1];
				node->quadMom[2] += node->children[i]->quadMom[2];
				node->quadMom[3] += node->children[i]->quadMom[3];
				node->quadMom[4] += node->children[i]->quadMom[4];
				node->quadMom[5] += node->children[i]->quadMom[5];
			}
		}
	}
}


/**
 * Given an (sub-)octree with root at node, insert
 * the body at position r into the octree, splitting nodes as needed.
 *
 * @param node, the root of the octree in which to insert the point.
 * @param r, the position of the body to insert.
 * @param m, the mass of the body to insert.
 */
void insertBodyOctree_NB(NBOctreeNode_t* node, const double* r, double m) {
	if (node == NULL) {
		return;
	}

	if (node->N == 0) {
		node->mass = m;
		node->com[0] = r[0];
		node->com[1] = r[1];
		node->com[2] = r[2];
	} else if (node->N > 1) {
		NBOctant_t oct = determineOctant_NB(node->center, r);
		if (node->children[oct]) {
			insertBodyOctree_NB(node->children[oct], r, m);
		} else {
			createNBodyChildNode(node, oct, r, m);
		}
	} else { /* N == 1 */
		pushNodePointToChild_NB(node);

		NBOctant_t oct = determineOctant_NB(node->center, r);
		if (node->children[oct]) {
			//if both node's point and the new point fall into same octant
			insertBodyOctree_NB(node->children[oct], r, m);
		} else {
			createNBodyChildNode(node, oct, r, m);
		}
	}

	//finally, increment before returning.
	++(node->N);
}

// void insertNeighbourBody_NB(NBOctreeNode_t* node, const double* r, double m) {
// 	if (node == NULL) {
// 		return;
// 	}

// 	if (node->N == 0) {
// 		node->mass = m;
// 		node->com[0] = r[0];
// 		node->com[1] = r[1];
// 		node->com[2] = r[2];
// 	} else if (node->N > 1) {
// 		NBOctant_t oct = determineOctant_NB(node->center, r);
// 		if (node->children[oct]) {
// 			//recurse
// 			insertNeighbourBody_NB(node->children[oct], r, m);
// 		} else {
// 			NBOctreeNode_t* newChild = createNBodyChildNode(node, oct, r, m);
// 			newChild->isNeighbourNode = 1;
// 		}
// 	} else { /* N == 1 */
// 		pushNodePointToChild_NB(node);

// 		NBOctant_t oct = determineOctant_NB(node->center, r);
// 		if (node->children[oct]) {
// 			//if both node's point and the new point fall into same octant, recurse
// 			insertNeighbourBody_NB(node->children[oct], r, m);
// 		} else {
// 			NBOctreeNode_t* newChild = createNBodyChildNode(node, oct, r, m);
// 			newChild->isNeighbourNode = 1;
// 		}
// 	}

// 	//finally, increment before returning.
// 	++(node->N);
// 	node->isNeighbourNode = 1;
// }



/**
 * Given a (sub-)Octree, free all nodes which are empty.
 *
 * @node, the root node of the (sub-)tree ot free.
 */
void pruneEmptyChildren_NB(NBOctreeNode_t* node) {
	if (node == NULL || node->children == NULL) {
		return;
	}

	int nullCount = 0;;
	for (int i = 0; i < 8; ++i) {
		if (node->children[i] == NULL) {
			++nullCount;
		} else if (node->children[i]->N == 0) {
			freeOctreeNode_NB(node->children[i]);
			node->children[i] = NULL;
			++nullCount;
		} else {
			pruneEmptyChildren_NB(node->children[i]);
		}
	}

	if (nullCount == 8) {
		free(node->children);
		node->children = NULL;
	}
}


/**
 * Given a list of N positions and masses, bounded by (-domainSize, domainSize)
 * in each dimension, create an octree containing all bodies.
 *
 * @param N, the number of bodies.
 * @param r, an array of 3*N values representing the bodies' positions.
 * @param m, an array of N values representing the bodies' masses.
 * @param domainSize, the size of the entire octree.
 * @return the newly created octree.
 */
NBOctree_t* buildOctree_NB(long N, const double* r, const double* m, double domainSize) {
	if (N == 0) {
		return NULL;
	}

	NBOctreeNode_t* root = createNBodyRootNode(domainSize);
	long i;
	for (i = 0; i < N; ++i) {
		insertBodyOctree_NB(root, r+(3*i), m[i]);
	}

	computeMassVals_NB(root);

	NBOctree_t* tree = (NBOctree_t*) malloc(sizeof(NBOctree_t));
	tree->root = root;
	return tree;
}


/**
 * Build an octree in place if possible.
 * Given a list of N positions and masses, bounded by (-domainSize, domainSize)
 * in each dimension, empty the octree pointed to tree_ptr and then
 * fill it with the N bodies.
 * If tree_ptr points to NULL, a new octree is created.
 * If tree_ptr points to a tree whose size is not domainSize, the existing tree
 * is freed and a new tree is created and returned in tree_ptr.
 *
 * @param N, the number of bodies.
 * @param r, an array of 3*N values representing the bodies' positions.
 * @param m, an array of N values representing the bodies' masses.
 * @param domainSize, the size of the entire octree.
 * @param tree_ptr, a pointer to the existing tree to update.
 * @return 1 iff the build actually occurred in place.
 */
int buildOctreeInPlace_NB(long N, const double* r, const double* m, double domainSize, NBOctree_t** tree_ptr) {
	if (tree_ptr == NULL) {
		return 0;
	}

	if (*tree_ptr == NULL) {
		*tree_ptr = buildOctree_NB(N, r, m, domainSize);
		return 0;
	}

	//TOOD fix below and actually do it in place.

//	NBOctreeNode_t* root = (*tree_ptr)->root;
//	if (fabs(domainSize - root->size) < _EPS) {
//		//domainSize and octree size are equal
//		emptyOctree_NB(root);
//	 	long i;
//	 	for (i = 0; i < N; ++i) {
//	 		insertBodyOctree_NB(root, r+(3*i), m[i]);
//	 	}
//	 	pruneEmptyChildren_NB(root);
//
//	 	computeMassVals_NB(root);
//	 	return 1;
//	} else {
		freeOctree_NB(*tree_ptr);
		*tree_ptr = buildOctree_NB(N, r, m, domainSize);
		return 0;
//	}
}


/**
 * Both dst and src should have children already, in order to maintain
 * linkage from parent node to child node without adding reverse links.
 */
void mergeOctreeNodes_NB(NBOctreeNode_t* dst, NBOctreeNode_t* src) {
	if (dst == NULL || src == NULL || dst->children == NULL) {
		return;
	}

	//assert(dst->center[0] == src->center[0]);
	//assert(dst->center[1] == src->center[1]);
	//assert(dst->center[2] == src->center[2]);
	//assert(dst->size = src->size);

	//4 cases to consider. src is always non-empty.
	// 1/ dst child is empty; src child is non-empty
	// 2/ dst child is leaf, src child is non-empty: insert dst leaf to src subtree and swap
	// 3/ dst child is internal; src child is internal: recurse the merge
	// 4/ dst child is internal; src child is leaf:

	dst->N = 0;
	for (int i = 0; i < 8; ++i) {
		if (dst->children[i] != NULL && src->children[i] != NULL) {
			if (dst->children[i]->children == NULL) {
				//case 2
				insertBodyOctree_NB(
					src->children[i],
					(dst->children[i]->com),
					dst->children[i]->mass
				);
				freeOctreeNode_NB(dst->children[i]);
				dst->children[i] = src->children[i];
				src->children[i] = NULL;
			} else if (src->children[i]->children == NULL) {
				//case 4;
				insertBodyOctree_NB(
					dst->children[i],
					(src->children[i]->com),
					src->children[i]->mass
				);
				//src freed at end so whatever with this node
			} else {
				//case 3
				mergeOctreeNodes_NB(dst->children[i], src->children[i]);
			}
			dst->N += dst->children[i]->N;
		} else if (src->children[i] != NULL) {
			//case 1;
			dst->children[i] = src->children[i];
			src->children[i] = NULL;
			dst->N += dst->children[i]->N;
			setOctChild_NB(dst->childByte, i);
		} else if (dst->children[i] != NULL) {
			dst->N += dst->children[i]->N;
		}
	}
}


/**
 * Merge two octrees in place, moving and inserting from srcTree to dstTree.
 * srcTree is left with many or all of its nodes removed and moved to dstTree.
 *
 * @param dstTree, the destination tree for the merge result
 * @param srcTree, the tree to merge cells from.
 */
void mergeOctreesInPlace_NB(NBOctree_t* dstTree, NBOctree_t* srcTree) {
	if (dstTree == NULL || srcTree == NULL || srcTree->root == NULL) {
		return;
	}

	if (dstTree->root == NULL) {
		dstTree->root = srcTree->root;
		return;
	}

	//dst is a single node; insert to src and swap
	if (dstTree->root->children == NULL) {
		insertBodyOctree_NB(srcTree->root, (dstTree->root->com), dstTree->root->mass);
		dstTree->root = srcTree->root;
		srcTree->root = NULL;
		return;
	}

	//src is single node; insert it to dst.
	if (srcTree->root->children == NULL) {
		insertBodyOctree_NB(dstTree->root, (srcTree->root->com), srcTree->root->mass);
		return;
	}

	//otherwise now both roots have children so we can enter the main algorithm
	mergeOctreeNodes_NB(dstTree->root, srcTree->root);
}


/**
 * Walk the octree start at root, determining the interactions
 * that should be including for the test point as position r.
 *
 * @param root, the roto of the tree to traverse.
 * @param r, an array of size 3 representing the position of the test point.
 * @param mac, a pointer to the multipole acceptance function.
 * @param theta, a parameter to the MAC.
 * @param walkList, an array of size at least root->N used as working space.
 * @param[out] interactList, an array of size at least root->N in which to
 *                           store the interacting cells.
 * @return the length of the interactList.
 */
long traverseTreeInteractionList_NB(const NBOctree_t* tree,
	const double* r,
	MAC_ptr mac,
	double theta,
	const NBOctreeNode_t** walkList,
	const NBOctreeNode_t** interactList)
{
	if (tree == NULL) {
		return 0;
	}

	const NBOctreeNode_t* node;
	walkList[0] = tree->root;
	long walkIdx = 1;
	long intIdx = 0;
	while (walkIdx > 0) {
		node = walkList[--walkIdx];
		if (node->children != NULL) {
			for (int i = 0; i < 8; ++i) {
				if (node->children[i] == NULL) {
					continue;
				}
				if ((*mac)(node->children[i], r, theta)) {
					interactList[intIdx++] = node->children[i];
				} else {
					walkList[walkIdx++] = node->children[i];
				}
			}
		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
			interactList[intIdx++] = node;
		}
	}

	return intIdx;
}


/**
 * Traverse a node recursively, building up the interaction list based on mac.
 */
void traverseNode(const NBOctreeNode_t* node, const double* r,
	MAC_ptr mac, double theta, const NBOctreeNode_t** interactList, long* intIdx) {
	if (node->children != NULL) {
		for(int i = 0; i < 8; ++i) {
			//this complex condition ensures the call containing r itself is not added.
			if(node->children[i] == NULL ||
				(node->children[i]->N == 1 && node->children[i]->com[0]==r[0] &&
				node->children[i]->com[1]==r[1] && node->children[i]->com[2]==r[2]))
			{
				continue;
			}
			if ((*mac)(node->children[i], r, theta)) {
				interactList[*intIdx] = node->children[i];
				++(*intIdx);
			} else {
				traverseNode(node->children[i], r, mac, theta, interactList, intIdx);
			}
		}
	} else {
		interactList[*intIdx] = node;
		++(*intIdx);
	}
}


/**
 * Walk the octree start at root, determining the interactions
 * that should be including for the test point as position r.
 * This method works recursively.
 *
 * @param root, the roto of the tree to traverse.
 * @param r, an array of size 3 representing the position of the test point.
 * @param mac, a pointer to the multipole acceptance function.
 * @param theta, a parameter to the MAC.
 * @param walkList, an array of size at least root->N used as working space.
 * @param[out] interactList, an array of size at least root->N in which to
 *                           store the interacting cells.
 * @return the length of the interactList.
 */
long traverseTreeInteractionListRec_NB(const NBOctree_t* tree,
	const double* r,
	MAC_ptr mac,
	double theta,
	const NBOctreeNode_t** walkList,
	const NBOctreeNode_t** interactList)
{
	if (tree == NULL) {
		return 0;
	}

	long intIdx = 0;
	traverseNode(tree->root, r, mac, theta, interactList, &intIdx);
	return intIdx;
}





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
	size_t* treeSize)
{
	if (tree == NULL || treeData == NULL || treeSize == NULL) {
		return;
	}

	//map
    buildOctreeInPlace_NB(localN, r, m, domainSize, tree);

    int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int ntrees = world_size;

    //reduce

	//get steps by bit-hacked ceil(log_2(ntrees))
	int nsteps = 0;
	int ntreetmp = ntrees;
	while (ntreetmp >>= 1) { ++nsteps; }
	if ((ntrees & (ntrees-1))) { ++nsteps; } //round up if not exactly a power of 2

	NBOctree_t* tmpTree = NULL;
	size_t newTreeSize;
	MPI_Status stat;
	int stepSize = 1;
	for (int i = 0; i < nsteps; ++i) {
		//rank k+stepsize sends to rank k to do the merge.
		for (int k = 0; k < ntrees; k += 2*stepSize) {
			if (world_rank == k && world_rank + stepSize < world_size) {
				// fprintf(stderr, "RANK %d RECEIVING FROM %d\n", world_rank, k + stepSize);
				MPI_Recv(&newTreeSize, 1, MPI_LONG,
					     k + stepSize, 0, MPI_COMM_WORLD, &stat);
				if (newTreeSize > *treeSize) {
					*treeData = realloc(*treeData, sizeof(uint8_t)*newTreeSize);
					*treeSize = newTreeSize;
				}
				MPI_Recv(*treeData, *treeSize, MPI_UINT8_T,
						 k + stepSize, 0, MPI_COMM_WORLD, &stat);

	            tmpTree = deserializeOctree_NB(*treeData, *treeSize);
	            mergeOctreesInPlace_NB(*tree, tmpTree);
				freeOctree_NB(tmpTree);

			} else if (world_rank == k + stepSize) {
				// fprintf(stderr, "RANK %d SENDING TO %d\n", world_rank, k);
				newTreeSize = serializeOctree_NB(*tree, treeData, treeSize);
				MPI_Send(&newTreeSize, 1, MPI_LONG,
						 k, 0, MPI_COMM_WORLD);
				MPI_Send(*treeData, newTreeSize, MPI_UINT8_T,
						 k, 0, MPI_COMM_WORLD);
			}
		}
		stepSize <<= 1;
	}

	if (world_rank == 0) {
		computeMassVals_NB((*tree)->root);
	}
}
#endif