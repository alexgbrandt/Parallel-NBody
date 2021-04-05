
#include "NBodySerialize.h"

#include <string.h>

size_t serializeHOTNodeList_NB(NBodyHOTNode_t* nodeHead, uint8_t** data_p, size_t* curAlloc) {
	if (nodeHead == NULL || data_p == NULL || curAlloc == NULL) {
		return 0;
	}

	NBodyHOTNode_t* node = nodeHead;

	int nnodes = 0;
	while (node != NULL) {
		++nnodes;
		node = node->next;
	}

	if (*data_p == NULL || *curAlloc < nnodes*HOTNodeSerSize) {
		*data_p = realloc(*data_p, sizeof(uint8_t)*nnodes*HOTNodeSerSize);
		*curAlloc = nnodes*HOTNodeSerSize;
	}

	uint8_t* data = *data_p;
	node = nodeHead;

	while (node != NULL) {
		memcpy(data, (void*) node, sizeof(uint8_t)*HOTNodeSerSize);
		node = node->next;
		data += HOTNodeSerSize;
	}

	return nnodes*HOTNodeSerSize;
}

size_t serializeHOTChildrenToList_NB(NBodyHOTNode_t* children[8], uint8_t** data_p, size_t* curAlloc) {
	if (data_p == NULL || curAlloc == NULL) {
		return 0;
	}

	int count = 0;
	for (int i = 0; i < 8; ++i) {
		count += (children[i] != NULL);
	}

	if (*data_p == NULL || *curAlloc < count*HOTNodeSerSize) {
		*data_p = realloc(*data_p, sizeof(uint8_t)*count*HOTNodeSerSize);
		*curAlloc = count*HOTNodeSerSize;
	}

	uint8_t* data = *data_p;
	for (int i = 0; i < 8; ++i) {
		if (children[i] != NULL) {
			memcpy(data, (void*) children[i], sizeof(uint8_t)*HOTNodeSerSize);
			data += HOTNodeSerSize;
		}
	}

	return count*HOTNodeSerSize;
}


NBodyHOTNode_t* deserializeHOTNodeList_NB(uint8_t* data, size_t size) {
	if (data == NULL || size == 0) {
		return NULL;
	}

	NBodyHOTNode_t* head = createHOTNode_NB();
	memcpy((void*) head, data, HOTNodeSerSize);
	NBodyHOTNode_t* tail = head;

	uint8_t* dataIter = data + HOTNodeSerSize;
	NBodyHOTNode_t* node;
	while (dataIter < data + size) {
		node = createHOTNode_NB();
		memcpy((void*) node, dataIter, sizeof(uint8_t)*HOTNodeSerSize);
		tail->next = node;
		tail = node;
		dataIter += HOTNodeSerSize;
	}

	return head;
}


/**
 * Serialize an array of 8 hashed octree node children.
 * Some entries in the array may be NULL, but they
 * are still serialized into null-data.
 * The resulting serialization is exactly HOTChildrenSerSize bytes.
 *
 * @param children: the array of HOT node children
 * @param data[out]: an array of bytes of size at least HOTChildrenSerSize
 */
void serializeHOTChildren_NB(NBodyHOTNode_t* children[8], uint8_t* data) {

	for (int i = 0; i < 8; ++i) {
		if (children[i] == NULL) {
			memset(data, 0, HOTNodeSerSize);
		} else {

			memcpy(data, children[i], HOTNodeSerSize);
		}
		NBodyHOTNode_t* node = (NBodyHOTNode_t*) data;
		// fprintf(stderr, "[ser] %d (%.5f, %.5f, %.5f), %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", i, node->center[0], node->center[1], node->center[2], node->key, node->N, node->com[0], node->com[1], node->com[2], node->mass, node->owningRank);
		data += HOTNodeSerSize;
	}
}

/**
 * Deserialize an array of 8 hashed octree node children.
 * Some resulting entries in the array may be NULL.
 * The data to deserialize is exactly HOTChildrenSerSize bytes.
 *
 * @param children[out]: the array of HOT node children
 * @param data: the array of bytes encoding the octree children.
 */
void deserializeHOTChildren_NB(NBodyHOTNode_t* children[8], uint8_t* data) {

	NBodyHOTNode_t* dst = createHOTNode_NB();
	for (int i = 0; i < 8; ++i) {
		children[i] = NULL;
		memcpy((void*) dst, data, HOTNodeSerSize);
		NBodyHOTNode_t* node = dst;
		// fprintf(stderr, "[des] %d (%.5f, %.5f, %.5f), %llo, N = %ld, COM=(%.5f, %.5f, %.5f), M=%.5f, RANK=%d\n", i, node->center[0], node->center[1], node->center[2], node->key, node->N, node->com[0], node->com[1], node->com[2], node->mass, node->owningRank);

		if (dst->N != 0) {
			children[i] = dst;
			dst = createHOTNode_NB();
		}
		data += HOTNodeSerSize;
	}
	freeHOTNode_NB(dst);
}


//pre-order traversal for easier de-serialization
uint8_t* _serializeOctNode_NB(NBOctreeNode_t* node, uint8_t* data) {
	size_t datachunk = sizeof(NBOctreeNode_t) - sizeof(NBOctreeNode_t**);

	//assumes that the list to children is at end of struct.
	//serailize everything else into data.
	memcpy(data, (void*) node, datachunk);
	data += datachunk;

	if (node->children != NULL) {
		for (int i = 0; i < 8; ++i) {
			if (node->children[i] != NULL) {
				data = _serializeOctNode_NB(node->children[i], data);
			}
		}
	}

	return data;
}


size_t serializeOctree_NB(NBOctree_t* tree, uint8_t** data_p, size_t* alloc) {
	if (data_p == NULL || tree == NULL) {
		return 0;
	}

	long nnodes = countNodes_NB(tree);

	size_t datachunk = sizeof(NBOctreeNode_t) - sizeof(NBOctreeNode_t**);

	uint8_t* datahead = *data_p;
	if (datahead == NULL) {
		datahead = (uint8_t*) malloc(datachunk*nnodes);
	} else if (*alloc < nnodes*datachunk) {
		datahead = (uint8_t*) realloc(datahead, datachunk*nnodes);
		*alloc = datachunk*nnodes;
	}

	_serializeOctNode_NB(tree->root, datahead);

	*data_p = datahead;
	return nnodes*datachunk;
}


uint8_t* _deserializeOctNode_NB(uint8_t* data, NBOctreeNode_t** retNode) {
	size_t datachunk = sizeof(NBOctreeNode_t) - sizeof(NBOctreeNode_t**);

	NBOctreeNode_t* node = (NBOctreeNode_t*) calloc(1, sizeof(NBOctreeNode_t));
	memcpy( (void*) node, data, datachunk);

	*retNode = node;
	data += datachunk;

	if (node->childByte == 0) {
		return data;
	}

	node->children = calloc(8, sizeof(NBOctreeNode_t*));
	for (int i = 0; i < 8; ++i) {
		if (hasOctChild_NB(node->childByte, i)) {
			data = _deserializeOctNode_NB(data, node->children + i);
		}
	}

	return data;
}


NBOctree_t* deserializeOctree_NB(uint8_t* data, size_t size) {
	if (data == NULL || size == 0) {
		return NULL;
	}

	NBOctree_t* tree = (NBOctree_t*) malloc(sizeof(NBOctree_t));

	_deserializeOctNode_NB(data, &(tree->root));

	return tree;
}


