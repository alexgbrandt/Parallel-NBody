
#include "NBodyHashedOctree.h"

#if NBODY_MPI
#include <mpi.h>

#include "NBodySerialize.h"
#endif

typedef enum HOTTraverseTags {
	KeyRequestTag = 0x17,
	ChildSizeTag,
	ChildDataTag
} HOTTraverseTags_t;



/**
 * Create an empty HOT node with default values.
 *
 * @return the newly constructed node.
 */
NBodyHOTNode_t* createHOTNode_NB() {
	NBodyHOTNode_t* node = (NBodyHOTNode_t*) malloc(sizeof(NBodyHOTNode_t));
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

	node->key = 0;
	node->childByte = 0;
	node->childrenLocal = 1;
	node->requested = 0;
	node->owningRank = 0;
	node->next = NULL;

	return node;
}

NBodyHOTNode_t* createHOTChildNode_NB(
		NBodyHOTNode_t* parent,
		NBOctant_t oct,
		double halfSize,
		const double* r,
		double m)
{
	NBodyHOTNode_t* childNode = createHOTNode_NB();
	childNode->size = halfSize;

	childNode->key = getOctChildKey_NB(parent->key, oct);

	childNode->center[0] = parent->center[0] + (halfSize*OctantDir[3*oct]);
	childNode->center[1] = parent->center[1] + (halfSize*OctantDir[3*oct+1]);
	childNode->center[2] = parent->center[2] + (halfSize*OctantDir[3*oct+2]);
	childNode->com[0] = r[0];
	childNode->com[1] = r[1];
	childNode->com[2] = r[2];
	childNode->mass = m;
	childNode->N = 1;
	return childNode;
}

NBodyHOTNode_t* createHotParentNode_NB(NBodyHOTNode_t* childNode) {
	NBodyHOTNode_t* parentNode = createHOTNode_NB();

	parentNode->size = childNode->size * 2;
	parentNode->key = getOctParentKey_NB(childNode->key);

	NBOctant_t oct = getOctantFromKey_NB(childNode->key);
	double halfSize = childNode->size;
	parentNode->center[0] = childNode->center[0] - (halfSize*OctantDir[3*oct]);
	parentNode->center[1] = childNode->center[1] - (halfSize*OctantDir[3*oct+1]);
	parentNode->center[2] = childNode->center[2] - (halfSize*OctantDir[3*oct+2]);
	parentNode->mass = childNode->mass;
	parentNode->N = childNode->N;
	return parentNode;
}

NBodyHOTNode_t* copyHOTNode_NB(NBodyHOTNode_t* node) {
	NBodyHOTNode_t* ret = (NBodyHOTNode_t*) malloc(sizeof(NBodyHOTNode_t));
	memcpy(ret, node, sizeof(NBodyHOTNode_t));

	ret->next = NULL; //don't keep references to existing hash table
	return ret;
}

void _insertBodyHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, const double* r, double m, int owning_rank) {
	if (hot.data == NULL) {
		return;
	}

	spatialKey_t curKey;
	if (node == NULL) {
		curKey = ROOT_KEY;
		node = getHOTNode_NB(hot, curKey);
	} else {
		curKey = node->key;
	}

	NBOctant_t oct;
	while (node != NULL) {
		if (node->N == 0) {
			//this case should only happen for the first insert to the root
			node->mass = m;
			node->com[0] = r[0];
			node->com[1] = r[1];
			node->com[2] = r[2];
			node->N = 1;
			node = NULL;
		} else if (node->N > 1) {
			oct = determineOctant_NB(node->center, r);
			node->N += 1;
			if (hasOctChild_NB(node->childByte, oct)) {
				curKey = getOctChildKey_NB(curKey, oct);
				node = getHOTNode_NB(hot, curKey);
			} else {
				NBodyHOTNode_t* child = createHOTChildNode_NB(node, oct, node->size*0.5, r, m);
				child->owningRank = owning_rank;
				setOctChild_NB(node->childByte, oct);
				insertHOTNode_NB(hot, child);
				node = NULL;
			}
		} else { /* N == 1 */
			pushHOTNodePointToChild_NB(hot, node);
			oct = determineOctant_NB(node->center, r);
			node->N += 1;
			if (hasOctChild_NB(node->childByte, oct)) {
				curKey = getOctChildKey_NB(curKey, oct);
				node = getHOTNode_NB(hot, curKey);
			} else {
				NBodyHOTNode_t* child = createHOTChildNode_NB(node, oct, node->size*0.5, r, m);
				child->owningRank = owning_rank;
				setOctChild_NB(node->childByte, oct);
				insertHOTNode_NB(hot, child);
				node = NULL;
			}
		}
	}
}


void insertBodyHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, const double* r, double m) {
	int owning_rank = 0;
#if NBODY_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &owning_rank);
#endif
    _insertBodyHOT_NB(hot, node, r, m, owning_rank);
}


NBodyHOT_t buildHashedOctreeInPlace_NB(long N, const double* r, const double* m, double domainSize, NBodyHOT_t hot)	{

	//TODO make this smarter maybe
	if (hot.data == NULL) {
		hot = createHashedOctree_NB(DEFAULT_HASHED_OCTREE_H);
	} else {
		clearHashedOctree_NB(hot);
	}

	int owning_rank = 0;
#if NBODY_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &owning_rank);
#endif

	NBodyHOTNode_t* root = createHOTNode_NB();
	root->size = domainSize;
	root->key = ROOT_KEY;
	root->owningRank = owning_rank;
	insertHOTNode_NB(hot, root);


	long i;
	for (i = 0; i < N; ++i) {
		_insertBodyHOT_NB(hot, NULL, r+(3*i), m[i], owning_rank);
	}

	if (N > 0) {
		computeMassValsHOT_NB(hot);
	}

	return hot;
}


void _computeMassValsHOTNode_NB(NBodyHOT_t hot, NBodyHOTNode_t* node) {
	//assert(node != NULL);
	//assert(node->N > 1);

	//reset current node before accumulating
	node->mass = 0.0;
	node->com[0] = 0.0;
	node->com[1] = 0.0;
	node->com[2] = 0.0;
	node->N = 0;

	//compute c.o.m. first.
	NBodyHOTNode_t* childNode;
	spatialKey_t childKey;
	double mass;
	int i;
	for (i = 0; i < 8; ++i) {
		if (hasOctChild_NB(node->childByte, i)) {
			childKey = getOctChildKey_NB(node->key, i);
			childNode = getHOTNode_NB(hot, childKey);
			mass = childNode->mass;
			node->mass += mass;
			node->N += childNode->N;
			node->com[0] += mass*childNode->com[0];
			node->com[1] += mass*childNode->com[1];
			node->com[2] += mass*childNode->com[2];

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
		if (hasOctChild_NB(node->childByte, i)) {
			childKey = getOctChildKey_NB(node->key, i);
			childNode = getHOTNode_NB(hot, childKey);

			//first, get child's pos in c.o.m. reference frame
			pos[0] = childNode->com[0] - node->com[0];
			pos[1] = childNode->com[1] - node->com[1];
			pos[2] = childNode->com[2] - node->com[2];

			mass = childNode->mass;
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
			if (childNode->N > 1) {
				node->quadMom[0] += childNode->quadMom[0];
				node->quadMom[1] += childNode->quadMom[1];
				node->quadMom[2] += childNode->quadMom[2];
				node->quadMom[3] += childNode->quadMom[3];
				node->quadMom[4] += childNode->quadMom[4];
				node->quadMom[5] += childNode->quadMom[5];
			}
		}
	}

}

/**
 * Given a node, inNode, in a hashed octree, hot,
 * updage all the mass and multipole values of its parent,
 * grandparent, etc.
 *
 * @param hot: the hashed octree
 * @param inNode: the cell from which to travel upwards and update nodes
 */
void _updateMassValsHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* inNode) {
	if (inNode == NULL || hot.data == NULL) {
		return;
	}


	spatialKey_t parentKey = getOctParentKey_NB(inNode->key);
	NBodyHOTNode_t* parentNode = getHOTNode_NB(hot, parentKey);
	while (parentNode != NULL) {
		_computeMassValsHOTNode_NB(hot, parentNode);
		if (parentKey == ROOT_KEY) {
			break;
		}

		parentKey = getOctParentKey_NB(parentKey);
		parentNode = getHOTNode_NB(hot, parentKey);
	}
}


void _computeMassValsHOT_NB(NBodyHOT_t hot, NBodyHOTNode_t* node) {
	if (node == NULL || hot.data == NULL || node->N == 1) {
		return;
	}

	//do a depth-first search
	//recursively update children first.
	int i;
	spatialKey_t childKey;
	for (i = 0; i < 8; ++i) {
		if (hasOctChild_NB(node->childByte, i)) {
			childKey = getOctChildKey_NB(node->key, i);
			_computeMassValsHOT_NB(hot, getHOTNode_NB(hot, childKey));
		}
	}

	//now, actually update this node's values
	_computeMassValsHOTNode_NB(hot, node);
}


/**
 * Walk the octree starting at the root, determining the interactions
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
long traverseHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	const NBodyHOTNode_t** walkList,
	const NBodyHOTNode_t** interactList)
{
	if (hot.data == NULL) {
		return 0;
	}

	const NBodyHOTNode_t* node;
	const NBodyHOTNode_t* child;
	walkList[0] = getHOTNode_NB(hot, ROOT_KEY);
	long walkIdx = 1;
	long intIdx = 0;
	while (walkIdx > 0) {/**
 * Walk the octree starting at the root, determining the interactions
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
long traverseHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	const NBodyHOTNode_t** walkList,
	const NBodyHOTNode_t** interactList)
{
	if (hot.data == NULL) {
		return 0;
	}

	const NBodyHOTNode_t* node;
	const NBodyHOTNode_t* child;
	walkList[0] = getHOTNode_NB(hot, ROOT_KEY);
	long walkIdx = 1;
	long intIdx = 0;
	while (walkIdx > 0) {
		node = walkList[--walkIdx];
		if (node->childByte != 0) {
			for (int i = 0; i < 8; ++i) {
				if (!hasOctChild_NB(node->childByte, i)) {
					continue;
				}
				child = getHOTNode_NB(hot, getOctChildKey_NB(node->key, i));
				if ((*mac)(child, r, theta)) {
					interactList[intIdx++] = child;
				} else {
					walkList[walkIdx++] = child;
				}
			}
		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
			interactList[intIdx++] = node;
		}
	}

	return intIdx;
}
		node = walkList[--walkIdx];
		if (node->childByte != 0) {
			for (int i = 0; i < 8; ++i) {
				if (!hasOctChild_NB(node->childByte, i)) {
					continue;
				}
				child = getHOTNode_NB(hot, getOctChildKey_NB(node->key, i));
				if ((*mac)(child, r, theta)) {
					interactList[intIdx++] = child;
				} else {
					walkList[walkIdx++] = child;
				}
			}
		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
			interactList[intIdx++] = node;
		}
	}

	return intIdx;
}


void moveHashedSubtree_NB(NBodyHOT_t dstHOT, NBodyHOT_t srcHOT, NBodyHOTNode_t* src) {
	if (src == NULL) {
		return;
	}

	NBodyHOTNode_t* srcChild;
	for (int i = 0; i < 8; ++i) {
		if (hasOctChild_NB(src->childByte, i)) {
			srcChild = getHOTNode_NB(srcHOT, getOctChildKey_NB(src->key, i));
			moveHashedSubtree_NB(dstHOT, srcHOT, srcChild);
		}
	}
	removeHOTNode_NB(srcHOT, src->key);
	insertHOTNode_NB(dstHOT, src);
}



/**
 * Both dst and src should have children already, in order to maintain
 * linkage from parent node to child node without adding reverse links.
 */
void _mergeHashedOctree_NB(NBodyHOT_t dstHOT, NBodyHOTNode_t* dst, NBodyHOT_t srcHOT, NBodyHOTNode_t* src) {
	if (dst == NULL || src == NULL || dst->childByte == 0) {
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
	NBodyHOTNode_t* dstChild;
	NBodyHOTNode_t* srcChild;
	for (int i = 0; i < 8; ++i) {
		if (hasOctChild_NB(dst->childByte, i) && hasOctChild_NB(src->childByte, i)) {
			dstChild = getHOTNode_NB(dstHOT, getOctChildKey_NB(dst->key, i));
			srcChild = getHOTNode_NB(srcHOT, getOctChildKey_NB(src->key, i));
			if (dstChild->childByte == 0) {
				//case 2
				insertBodyHOT_NB(srcHOT, srcChild, dstChild->com, dstChild->mass);

				removeHOTNode_NB(dstHOT, dstChild->key);
				freeHOTNode_NB(dstChild);

				moveHashedSubtree_NB(dstHOT, srcHOT, srcChild);
				unsetOctChild_NB(src->childByte, i);

				dstChild = getHOTNode_NB(dstHOT, getOctChildKey_NB(dst->key, i));
			} else if (srcChild->childByte == 0) {
				//case 4;
				insertBodyHOT_NB(dstHOT, dstChild, srcChild->com, srcChild->mass);
				removeHOTNode_NB(srcHOT, srcChild->key);
				unsetOctChild_NB(src->childByte, i);
				freeHOTNode_NB(srcChild);
			} else {
				//case 3
				_mergeHashedOctree_NB(dstHOT, dstChild, srcHOT, srcChild);
			}
			dst->N += dstChild->N;
		} else if (hasOctChild_NB(src->childByte, i)) {
			//case 1;
			//take srcChild and make it dstChild
			srcChild = getHOTNode_NB(srcHOT, getOctChildKey_NB(src->key, i));
			moveHashedSubtree_NB(dstHOT, srcHOT, srcChild);
			unsetOctChild_NB(src->childByte, i);

			dst->N += srcChild->N;
			setOctChild_NB(dst->childByte, i);
		} else if (hasOctChild_NB(dst->childByte, i)) {
			dstChild = getHOTNode_NB(dstHOT, getOctChildKey_NB(dst->key, i));
			dst->N += dstChild->N;
		}
	}
}

/**
 * Merge two octrees in place, moving and inserting from srcTree to dstTree.
 * srcTree is emptied and free'd by there operation and should not be used afterward.
 * dstTree and srcTree must have same hash table size.
 *
 * @param dstTree, a pointer to the destination tree for the merge result
 * @param srcTree, the tree to merge cells from.
 */
void mergeHashedOctrees_NB(NBodyHOT_t* dstTree, NBodyHOT_t srcTree) {
	if (dstTree == NULL || dstTree->data == NULL || srcTree.data == NULL || dstTree->h != srcTree.h) {
		return;
	}

	NBodyHOTNode_t* srcRoot;
	if ((srcRoot = getHOTNode_NB(srcTree, ROOT_KEY)) == NULL) {
		return;
	}

	NBodyHOTNode_t* dstRoot;
	if ((dstRoot = getHOTNode_NB(*dstTree, ROOT_KEY)) == NULL) {
		free(dstTree->data);
		dstTree->data = srcTree.data;
		return;
	}

	//dst is a single node; insert to src and swap
	if (dstRoot->childByte == 0) {
		insertBodyHOT_NB(srcTree, srcRoot, dstRoot->com, dstRoot->mass);
		freeHOTNode_NB(dstRoot);
		free(dstTree->data);
		dstTree->data = srcTree.data;
		return;
	}

	//src is single node; insert it to dst.
	if (srcRoot->childByte == 0) {
		insertBodyHOT_NB(*dstTree, dstRoot, srcRoot->com, srcRoot->mass);
		freeHOTNode_NB(srcRoot);
		free(srcTree.data);
		return;
	}

	//otherwise now both roots have children so we can enter the main algorithm
	_mergeHashedOctree_NB(*dstTree, dstRoot, srcTree, srcRoot);

	free(srcTree.data);
}


int _getBranchNodes_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, int world_rank, NBodyHOTNode_t** branchList_p) {
	if (node == NULL || hot.data == NULL || branchList_p == NULL) {
		return 0;
	}

	if (node->N == 1) {
		if (world_rank == node->owningRank) {
			*branchList_p = node;
			return 0;
		} else {
			return 1;
		}
	}

	//therefore N > 1 and node has children

	NBodyHOTNode_t* retHead = NULL;
	NBodyHOTNode_t* retTail = NULL;

	NBodyHOTNode_t* childBranches[8];
	int containsNeighbour = 0;
	spatialKey_t childKey;
	for (int i = 0; i < 8; ++i) {
		childBranches[i] = NULL;
		if (hasOctChild_NB(node->childByte, i)) {
			childKey = getOctChildKey_NB(node->key, i);
			int hasNeighbour = _getBranchNodes_NB(hot, getHOTNode_NB(hot, childKey), world_rank, childBranches + i);
			containsNeighbour |= hasNeighbour;

			if (hasNeighbour && childBranches[i] != NULL) {
				//then childBranches[i] is truly a list of childBranches
				if (retHead == NULL) {
					retHead = childBranches[i];
					retTail = childBranches[i];
				} else {
					retTail->next = childBranches[i];
				}
				while (childBranches[i]->next != NULL) {
					childBranches[i] = childBranches[i]->next;
				}
				retTail = childBranches[i];

				//use this to say we already added this child's branch nodes to the LL
				childBranches[i] = NULL;
			}
		}
	}

	if (containsNeighbour) {
		//in this case, we have a child with the neighbour,
		//and thus some branch nodes, and possibly some other children
		//which are not yet labelled as branch nodes.
		//The true branch nodes for the child which contains the neighbour
		//has already been added to the retHead linked list by previous loop.
		//Here, we create branch nodes out of any other children and add it to the LL.
		//If we are in the special case that this node is the parent of the leaf
		//node containing the neighbour, then this loop also creates the retHead LL.
		for (int i = 0; i < 8; ++i) {
			if (childBranches[i] != NULL) {
				NBodyHOTNode_t* tmp = copyHOTNode_NB(childBranches[i]);
				if (retHead == NULL) {
					retHead = tmp;
					retTail = tmp;
				} else {
					//don't need to iterate over tail because we are adding
					//exactly one node
					retTail->next = tmp;
					retTail = tmp;
				}
			}
		}
		*branchList_p = retHead;
		return 1;
	} else {
		//all childBranches are not really branches and this node itself is a
		//candidate for a branch node
		*branchList_p = node;
		return 0;
	}
}


NBodyHOTNode_t* getBranchNodes_NB(NBodyHOT_t hot) {
	if (hot.data == NULL) {
		return NULL;
	}

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	NBodyHOTNode_t* root = getHOTNode_NB(hot, ROOT_KEY);
	NBodyHOTNode_t* branchList = NULL;

	_getBranchNodes_NB(hot, root, world_rank, &branchList);

	return branchList;
}

#if NBODY_MPI
void exchangeBorderBodies_NB(NBodyHOT_t hot, long localN, double* r, double* m) {

	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double leftData[4];
    double rightData[4];
    leftData[0] = r[0];
    leftData[1] = r[1];
    leftData[2] = r[2];
    leftData[3] = m[0];
    rightData[0] = r[3*localN - 3];
    rightData[1] = r[3*localN - 2];
    rightData[2] = r[3*localN - 1];
    rightData[3] = m[localN - 1];

    int firstExchangePartner, secondExchangePartner;
    int firstBuffRight;
    if (world_rank % 2 == 0) {
        //even procs send right first, then send left second.
        firstExchangePartner = world_rank + 1 < world_size ? world_rank + 1 : -1;
        secondExchangePartner = world_rank - 1 >= 0 ? world_rank - 1 : -1;
        firstBuffRight = 1;
    } else {
        firstExchangePartner = world_rank - 1 >= 0 ? world_rank - 1 : -1;
        secondExchangePartner = world_rank + 1 < world_size ? world_rank + 1 : -1;
        firstBuffRight = 0;
    }


    MPI_Status stat;
    if (firstExchangePartner >= 0) {
        if (firstBuffRight) {
            MPI_Sendrecv_replace(rightData, 4, MPI_DOUBLE, firstExchangePartner, 0,
            					 firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
        } else {
        	MPI_Sendrecv_replace(leftData, 4, MPI_DOUBLE, firstExchangePartner, 0,
        						 firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
        }
    }
    if (secondExchangePartner >= 0) {
    	if (firstBuffRight) {
            MPI_Sendrecv_replace(leftData, 4, MPI_DOUBLE, secondExchangePartner, 0,
            					 secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
    	} else {
            MPI_Sendrecv_replace(rightData, 4, MPI_DOUBLE, secondExchangePartner, 0,
                                 secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
    	}
    }

    if (firstExchangePartner >= 0) {
        if (firstBuffRight) {
			_insertBodyHOT_NB(hot, NULL, rightData, rightData[3], firstExchangePartner);
		} else {
			_insertBodyHOT_NB(hot, NULL, leftData, leftData[3], firstExchangePartner);
		}
    }
    if (secondExchangePartner >= 0) {
    	if (firstBuffRight) {
			_insertBodyHOT_NB(hot, NULL, leftData, leftData[3], secondExchangePartner);
    	} else {
			_insertBodyHOT_NB(hot, NULL, rightData, rightData[3], secondExchangePartner);
    	}
    }

}

/**
 * Given a branch node, node, insert it into the hashed octree, hot,
 * and create "fill" nodes so that a path exists from the root
 * of the tree down to that node.
 */
void insertAndFillHOTNode_NB(NBodyHOT_t hot, NBodyHOTNode_t* inNode) {
	if (inNode == NULL || hot.data == NULL) {
		return;
	}

	spatialKey_t childKey, parentKey;
	NBOctant_t oct;

	NBodyHOTNode_t* node = inNode;
	NBodyHOTNode_t* parentNode;
	while (node != NULL) {
		childKey = node->key;
		insertHOTNode_NB(hot, node);

		oct = getOctantFromKey_NB(childKey);
		parentKey = getOctParentKey_NB(childKey);
		parentNode = getHOTNode_NB(hot, parentKey);
		if (parentNode == NULL) {
			parentNode = createHotParentNode_NB(node);
			setOctChild_NB(parentNode->childByte, oct);
			node = parentNode;
		} else {
			setOctChild_NB(parentNode->childByte, oct);
			node = NULL;
		}
	}

	_updateMassValsHOT_NB(hot, inNode);

}

void exchangeBranchNodes_NB(NBodyHOT_t hot,	uint8_t** branchData, size_t* branchAlloc)
{
	if (hot.data == NULL || branchData == NULL || branchAlloc == NULL) {
		return;
	}


	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    NBodyHOTNode_t* branchNodes = getBranchNodes_NB(hot);

	// for (int k = 0; k < world_size; ++k) {
	// 	if (world_rank == k) {
	// 		fprintf(stderr, "\nRANK %d\n", world_rank);
	// 		printHashedOctree_NB(hot);
	// 	}
	// 	usleep(1000000);
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }


    //Every node needs to broadcast its branch nodes
    NBodyHOTNode_t* otherBranches;
    size_t recvSize;
    for (int i = 0; i < world_size; ++i) {
    	if (world_rank == i) {
			recvSize = serializeHOTNodeList_NB(branchNodes, branchData, branchAlloc);
    	}
    	MPI_Bcast(&recvSize, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_COMM_WORLD);

    	if (world_rank != i) {
    		if (recvSize > *branchAlloc) {
    			*branchAlloc = recvSize;
    			*branchData = (uint8_t*) realloc(*branchData, sizeof(uint8_t)*(*branchAlloc));
    		}

    	}
    	MPI_Bcast(*branchData, recvSize, MPI_UINT8_T, i, MPI_COMM_WORLD);


    	if (world_rank != i) {
			otherBranches = deserializeHOTNodeList_NB(*branchData, recvSize);

			NBodyHOTNode_t* node = otherBranches;
            NBodyHOTNode_t *tmp, *existingNode;
            while (node != NULL) {
            	existingNode = getHOTNode_NB(hot, node->key);
            	if (existingNode != NULL) {
            		if (existingNode->owningRank != i) {
            			fprintf(stderr, "HMMM.. branch nodes incorrect??\n");
    					fprintf(stderr, "\nRANK %d\n", world_rank);
						fprintf(stderr, "received\n");
                		fprintf(stderr, "(%.5f, %.5f, %.5f), %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", node->center[0], node->center[1], node->center[2], node->key, node->N, node->com[0], node->com[1], node->com[2], node->mass, node->owningRank);
						fprintf(stderr, "had\n");
                		fprintf(stderr, "(%.5f, %.5f, %.5f), %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", existingNode->center[0], existingNode->center[1], existingNode->center[2], existingNode->key, existingNode->N, existingNode->com[0], existingNode->com[1], existingNode->com[2], existingNode->mass, existingNode->owningRank);
            			exit(1);
            		} else {
            			removeHOTNode_NB(hot, node->key);
            		}
            	}

            	tmp = node->next;
            	node->next = NULL;
            	node->childrenLocal = 0;
            	node->requested = 0;
				insertAndFillHOTNode_NB(hot, node);
            	node = tmp;
            }
    	}
    }
}


void requestAndInsertHOTChildren_NB(NBodyHOT_t hot, NBodyHOTNode_t* node, uint8_t** treeData, size_t* treeAlloc) {
	if (node == NULL || node->childrenLocal == 1) {
		return;
	}

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int other_rank = node->owningRank;
	spatialKey_t reqKey = node->key;

	MPI_Request req;
	MPI_Issend(&reqKey, 1, MPI_UNSIGNED_LONG_LONG, other_rank, KeyRequestTag, MPI_COMM_WORLD, &req);
	int done = 0;
	while (!done) {
		//avoid deadlock by first checking for incoming requests
		fulfillHOTChildrenRequests_NB(hot, treeData, treeAlloc);

		MPI_Test(&req, &done, MPI_STATUS_IGNORE);
	}

	size_t recvCount;
	MPI_Status stat;
	MPI_Recv(&recvCount, 1, MPI_UNSIGNED_LONG_LONG, other_rank, ChildSizeTag, MPI_COMM_WORLD, &stat);

	if (recvCount == 0) {
		//request key doesn't exist on that rank... what's happening?
		return;
	}
	if (recvCount > *treeAlloc) {
		*treeAlloc = recvCount;
		*treeData = (uint8_t*) realloc(*treeData, sizeof(uint8_t)*(*treeAlloc));
	}

	MPI_Recv(*treeData, recvCount, MPI_UINT8_T, other_rank, ChildDataTag, MPI_COMM_WORLD, &stat);

	NBodyHOTNode_t* children = deserializeHOTNodeList_NB(*treeData, recvCount);
	NBodyHOTNode_t *tmp;

	while (children != NULL) {
		tmp = children->next;
		children->next = NULL;
		children->childrenLocal = 0;
		insertHOTNode_NB(hot, children);

		children = tmp;
	}
	node->childrenLocal = 1;
}


void asyncFulfillHOTChildrenRequests_NB(NBodyHOT_t hot, uint8_t* treeData) {

	int hasMessage = 0;
	MPI_Status stat;

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	MPI_Iprobe(MPI_ANY_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &hasMessage, &stat);

	spatialKey_t reqKey;
	NBodyHOTNode_t* parent;
	NBodyHOTNode_t* children[8];
	while (hasMessage) {
		//get requested key
		MPI_Recv(&reqKey, 1, MPI_UNSIGNED_LONG_LONG, stat.MPI_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &stat);

		//see if it actually exists
		parent = getHOTNode_NB(hot, reqKey);
		size_t sendCount = 0;
		if (parent == NULL) {
			memset(treeData, 0, HOTChildrenSerSize);
			MPI_Send(treeData, HOTChildrenSerSize, MPI_UINT8_T, stat.MPI_SOURCE, ChildDataTag, MPI_COMM_WORLD);
		} else {
			//serialize extant children
			for (int i = 0; i < 8; ++i) {
				children[i] = NULL;
				if (hasOctChild_NB(parent->childByte, i)) {
					children[i] = getHOTNode_NB(hot, getOctChildKey_NB(parent->key, i));
				}
			}
			serializeHOTChildren_NB(children, treeData);
			MPI_Send(treeData, HOTChildrenSerSize, MPI_UINT8_T, stat.MPI_SOURCE, ChildDataTag, MPI_COMM_WORLD);
		}

		//check for another
		hasMessage = 0;
		MPI_Iprobe(MPI_ANY_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &hasMessage, &stat);
	}
}


void asyncRequestHOTChildren_NB(NBodyHOT_t hot, NBodyDeferData_t* deferData) {
	//before we make a reqest, lets cooperate and see if we can fulfill any incoming ones
	asyncFulfillHOTChildrenRequests_NB(hot, deferData->recvBuffer);

	if (deferData == NULL || deferData->parentNode->childrenLocal == 1 || deferData->parentNode->requested == 1) {
		return;
	}

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int other_rank = deferData->parentNode->owningRank;
	spatialKey_t reqKey = deferData->parentNode->key;
	deferData->parentNode->requested = 1;

	MPI_Isend(&reqKey, 1, MPI_UNSIGNED_LONG_LONG, other_rank, KeyRequestTag, MPI_COMM_WORLD, deferData->reqs);
	MPI_Irecv(deferData->recvBuffer, HOTChildrenSerSize, MPI_UINT8_T, other_rank, ChildDataTag, MPI_COMM_WORLD, deferData->reqs + 1);

}

void asyncInsertHOTChildren_NB(NBodyHOT_t hot, NBodyDeferData_t* deferData, uint8_t* treeData) {

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// fprintf(stderr, "RANK %d waiting to insert\n", world_rank);

	int done = 0;
	while (!done) {
		//avoid deadlock by first checking for incoming requests
		asyncFulfillHOTChildrenRequests_NB(hot, treeData);
		MPI_Testall(2, deferData->reqs, &done, MPI_STATUS_IGNORE);
	}

	MPI_Waitall(2, deferData->reqs, MPI_STATUS_IGNORE);
	// fprintf(stderr, "RANK %d done waiting\n", world_rank );


	NBodyHOTNode_t* children[8];
	deserializeHOTChildren_NB(children, deferData->recvBuffer);

	for (int i = 0; i < 8; ++i) {
		if (children[i] != NULL) {
			children[i]->childrenLocal = 0;
			children[i]->requested = 0;
			insertHOTNode_NB(hot, children[i]);
		}
	}
	deferData->parentNode->childrenLocal = 1;
	deferData->parentNode->requested = 0;
}

void fulfillHOTChildrenRequests_NB(NBodyHOT_t hot, uint8_t** treeData, size_t* treeAlloc) {

	int hasMessage = 0;
	MPI_Status stat;

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	MPI_Iprobe(MPI_ANY_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &hasMessage, &stat);

	spatialKey_t reqKey;
	NBodyHOTNode_t* parent;
	NBodyHOTNode_t* children[8];
	while (hasMessage) {
		//get requested key
		MPI_Recv(&reqKey, 1, MPI_UNSIGNED_LONG_LONG, stat.MPI_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &stat);

		//see if it actually exists
		parent = getHOTNode_NB(hot, reqKey);
		size_t sendCount = 0;
		if (parent == NULL) {
			MPI_Send(&sendCount, 1, MPI_UNSIGNED_LONG_LONG, stat.MPI_SOURCE, ChildSizeTag, MPI_COMM_WORLD);
		} else {
			//serialize extant children
			for (int i = 0; i < 8; ++i) {
				children[i] = NULL;
				if (hasOctChild_NB(parent->childByte, i)) {
					children[i] = getHOTNode_NB(hot, getOctChildKey_NB(parent->key, i));
				}
			}
			sendCount = serializeHOTChildrenToList_NB(children, treeData, treeAlloc);

			//send num bytes first and then the data
			MPI_Send(&sendCount, 1, MPI_UNSIGNED_LONG_LONG, stat.MPI_SOURCE, ChildSizeTag, MPI_COMM_WORLD);
			MPI_Send(*treeData, sendCount, MPI_UINT8_T, stat.MPI_SOURCE, ChildDataTag, MPI_COMM_WORLD);
		}

		//check for another
		hasMessage = 0;
		MPI_Iprobe(MPI_ANY_SOURCE, KeyRequestTag, MPI_COMM_WORLD, &hasMessage, &stat);
	}
}


/**
 * Walk the octree starting at the root, determining the interactions
 * that should be including for the test point as position r.
 * The hashed octree is not necessarily fully present on the local processor,
 * requests for non-local data may be required.
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
long traverseDistributedHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	NBodyHOTNode_t** walkList,
	NBodyHOTNode_t** interactList,
	uint8_t** treeData,
	size_t* treeAlloc)
{
	if (hot.data == NULL) {
		return 0;
	}

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	NBodyHOTNode_t* node;
	NBodyHOTNode_t* child;
	walkList[0] = getHOTNode_NB(hot, ROOT_KEY);
	long walkIdx = 1;
	long intIdx = 0;
	while (walkIdx > 0) {
		node = walkList[--walkIdx];

		if (node->childByte != 0 && !node->childrenLocal) {
			requestAndInsertHOTChildren_NB(hot, node, treeData, treeAlloc);
		}

		if (node->childByte != 0) {
			for (int i = 0; i < 8; ++i) {
				if (!hasOctChild_NB(node->childByte, i)) {
					continue;
				}

				child = getHOTNode_NB(hot, getOctChildKey_NB(node->key, i));
				if ((*mac)(child, r, theta)) {
					interactList[intIdx++] = child;
				} else {
					walkList[walkIdx++] = child;
				}
			}
		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
			interactList[intIdx++] = node;
		}
	}

	return intIdx;
}


long traverseDistribAsyncHOTInteractionList_NB(NBodyHOT_t hot,
	const double* r,
	HOTMAC_ptr mac,
	double theta,
	NBodyHOTNode_t** walkList,
	// NBodyHOTNode_t** deferList,
	NBodyHOTNode_t** interactList,
	uint8_t** treeData,
	size_t* treeAlloc)
{
	if (hot.data == NULL) {
		return 0;
	}

	int deferListSize = 40;
	long deferHead = 0; //always points to an empty slot in the defer list.
	long deferTail = 0;
	NBodyDeferData_t deferList[deferListSize];

	if (*treeAlloc < deferListSize*HOTChildrenSerSize) {
		*treeData = (uint8_t*) realloc(*treeData, deferListSize*HOTChildrenSerSize);
	}

	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	NBodyHOTNode_t* node;
	NBodyHOTNode_t* child;
	walkList[0] = getHOTNode_NB(hot, ROOT_KEY);
	long walkIdx = 1;
	long intIdx = 0;
	while (walkIdx > 0 || deferHead != deferTail) {
		//try to take off the walk list, otherwise receive data from the oldest defer
		//and go with that.
		if (walkIdx > 0) {
			node = walkList[--walkIdx];
		} else {
			// fprintf(stderr, "RANK %d walkList 0, deferTail: %d\n", world_rank, deferTail);
			//use deferHead as temp data for sending buffer during async insert to fulfill while waiting
			asyncInsertHOTChildren_NB(hot, deferList + deferTail, (*treeData) + (deferHead*HOTChildrenSerSize));
			node = deferList[deferTail].parentNode;
			deferTail = (deferTail + 1) % deferListSize;
		}

		if (node->childByte != 0 && !node->childrenLocal) {
			//+1 to keep at least one empty slot for the send buffer
			if ( ((deferHead + 1) % deferListSize) == deferTail) {
				// fprintf(stderr, "RANK %d head matches tail, deferTail: %d\n", world_rank, deferTail);
				asyncInsertHOTChildren_NB(hot, deferList + deferTail, (*treeData) + (deferHead*HOTChildrenSerSize));
				walkList[walkIdx++] = deferList[deferTail].parentNode;
				deferTail = (deferTail + 1) % deferListSize;
			}

			// fprintf(stderr, "RANK %d inserting into deferHead: %d\n", world_rank, deferHead);
			deferList[deferHead].parentNode = node;
			deferList[deferHead].recvBuffer = (*treeData) + (deferHead*HOTChildrenSerSize);
			asyncRequestHOTChildren_NB(hot, deferList + deferHead);
			// fprintf(stderr, "RANK %d finished insert to deferHead: %d\n", world_rank, deferHead );
			deferHead = (deferHead + 1) % deferListSize;
			continue;
		}

		if (node->childByte != 0) {
			for (int i = 0; i < 8; ++i) {
				if (!hasOctChild_NB(node->childByte, i)) {
					continue;
				}

				child = getHOTNode_NB(hot, getOctChildKey_NB(node->key, i));
				// if (child == NULL) {
				// 	fprintf(stderr, "\n\n#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nRANK %d CHILD WAS NULL %d, %llo, %llo\n", world_rank, i, node->key,  getOctChildKey_NB(node->key, i));
				// 	exit(1);
				// }
				if ((*mac)(child, r, theta)) {
					interactList[intIdx++] = child;
				} else {
					walkList[walkIdx++] = child;
				}
			}
		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
			interactList[intIdx++] = node;
		}
	}

	//for good measure, fulfill requests before returning
	asyncFulfillHOTChildrenRequests_NB(hot, (*treeData) + (deferHead*HOTChildrenSerSize));

	// fprintf(stderr, "RANK %d DONE LOOP\n", world_rank);
	// int allDone = 0;
	// MPI_Request req;
	// MPI_Ibarrier(MPI_COMM_WORLD, &req);
	// while (!allDone) {
	// 	asyncFulfillHOTChildrenRequests_NB(hot, *treeData);
	// 	MPI_Test(&req, &allDone, MPI_STATUS_IGNORE);
	// }
	// fprintf(stderr, "RANK %d DONE BARRIER\n\n", world_rank);


	return intIdx;
}





// long traverseDistribAsyncHOTInteractionList_NB(NBodyHOT_t hot,
// 	const double* r,
// 	HOTMAC_ptr mac,
// 	double theta,
// 	NBodyHOTNode_t** walkList,
// 	// NBodyHOTNode_t** deferList,
// 	NBodyHOTNode_t** interactList,
// 	uint8_t** treeData,
// 	size_t* treeAlloc)
// {
// 	if (hot.data == NULL) {
// 		return 0;
// 	}

// 	int deferListSize = 40;
// 	long deferHead = 0;
// 	long deferTail = 0;
// 	NBodyDeferData_t* deferList = malloc(sizeof(NBodyDeferData_t)*deferListSize);
// 	if (*treeAlloc < deferListSize*HOTChildrenSerSize) {
// 		*treeData = (uint8_t*) realloc(*treeData, deferListSize*HOTChildrenSerSize);
// 	}


// 	int world_rank, world_size;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


// 	for (int k = 0; k < world_size; ++k) {
// 		if (world_rank == k) {
// 			fprintf(stderr, "\nRANK %d\n", world_rank);
// 			printHashedOctree_NB(hot);
// 		}
// 		usleep(1000000);
// 		MPI_Barrier(MPI_COMM_WORLD);
// 	}

// 	uint8_t sendBuff[HOTChildrenSerSize];

// 	NBodyHOTNode_t* node;
// 	NBodyHOTNode_t* child;
// 	walkList[0] = getHOTNode_NB(hot, ROOT_KEY);
// 	long walkIdx = 1;
// 	long intIdx = 0;
// 	while (walkIdx > 0 || deferHead != deferTail) {
// 		// fprintf(stderr, "RANK %d walkIdx %d, deferHead %d, deferTail %d\n", world_rank, walkIdx, deferHead, deferTail);
// 		//try to take off the walk list, otherwise receive data from the oldest defer
// 		//and go with that.
// 		if (walkIdx > 0) {
// 			node = walkList[--walkIdx];
// 		} else {
// 			// fprintf(stderr, "RANK %d walkList 0, deferTail: %d\n", world_rank, deferTail);
// 			//use deferHead as temp data for sending buffer during async insert to fulfill while waiting
// 			asyncInsertHOTChildren_NB(hot, deferList + deferTail, sendBuff);
// 			node = deferList[deferTail].parentNode;
// 			deferTail = (deferTail + 1) % deferListSize;
// 		}

// 		if (node->childByte != 0 && !node->childrenLocal) {
// 			if ( ((deferHead + 1) % deferListSize) == deferTail) {
// 				// fprintf(stderr, "RANK %d head matches tail, deferTail: %d\n", world_rank, deferTail);
// 				asyncInsertHOTChildren_NB(hot, deferList + deferTail, sendBuff);
// 				walkList[walkIdx++] = deferList[deferTail].parentNode;
// 				deferTail = (deferTail + 1) % deferListSize;
// 			}

// 			// fprintf(stderr, "RANK %d inserting into deferHead: %d\n", world_rank, deferHead);
// 			deferList[deferHead].parentNode = node;
// 			deferList[deferHead].recvBuffer = (*treeData) + (deferHead*HOTChildrenSerSize);
// 			asyncRequestHOTChildren_NB(hot, deferList + deferHead);
// 			// fprintf(stderr, "RANK %d finished insert to deferHead: %d\n", world_rank, deferHead );
// 			deferHead = (deferHead + 1) % deferListSize;
// 			continue;
// 		}

// 		if (node->childByte != 0) {
// 			for (int i = 0; i < 8; ++i) {
// 				if (!hasOctChild_NB(node->childByte, i)) {
// 					continue;
// 				}

// 				child = getHOTNode_NB(hot, getOctChildKey_NB(node->key, i));
// 				// if (child == NULL) {
// 				// 	fprintf(stderr, "\n\n#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nRANK %d CHILD WAS NULL %d, %llo, %llo\n", world_rank, i, node->key,  getOctChildKey_NB(node->key, i));
// 				// 	exit(1);
// 				// }
// 				if ((*mac)(child, r, theta)) {
// 					interactList[intIdx++] = child;
// 				} else {
// 					walkList[walkIdx++] = child;
// 				}
// 			}
// 		} else if (node->com[0] != r[0] || node->com[1] != r[1] || node->com[2] != r[2] ) {
// 			interactList[intIdx++] = node;
// 		}
// 	}

// 	//for good measure, fulfill requests before returning
// 	asyncFulfillHOTChildrenRequests_NB(hot, *treeData);

// 	// fprintf(stderr, "RANK %d DONE LOOP\n", world_rank);
// 	// MPI_Barrier(MPI_COMM_WORLD);
// 	// fprintf(stderr, "RANK %d DONE BARRIER\n\n", world_rank);


// 	return intIdx;
// }


#endif