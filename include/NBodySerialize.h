
#ifndef _NBODY_SERIALIZE
#define _NBODY_SERIALIZE_

#include <stdlib.h>
#include <stdint.h>

#include "NBodyOctree.h"
#include "NBodyHashedOctree.h"

#ifdef __cplusplus
extern "C" {
#endif

#define HOTNodeSerSize (sizeof(NBodyHOTNode_t) - sizeof(NBodyHOTNode_t*))
#define HOTChildrenSerSize (8*(HOTNodeSerSize))


/**
 * Serialize a linked list of hashed octree nodes, given by nodeHead.
 *
 * @param nodeHead: the head of the linked list to serialize
 * @param data_p[in,out]: a pointer to a possibly preallocated array of bytes holding the serialization
 * @param curAlloc[in,out]: a pointer to the number of bytes currently allocated in the data array;
 *                          returns in the pointer the new allocation if reallocated.
 * @return the number of bytes of the serialization
 */
size_t serializeHOTNodeList_NB(NBodyHOTNode_t* nodeHead, uint8_t** data_p, size_t* curAlloc);


/**
 * Serialize an array of 8 hashed octree node children.
 * Some entries in the array may be NULL, thus the returned size
 * may be less than 8 * the size of a HOT node.
 *
 * @param children: the array of HOT node children
 * @param data_p[in,out]: a pointer to a possibly preallocated array of bytes holding the serialization
 * @param curAlloc[in,out]: a pointer to the number of bytes currently allocated in the data array;
 *                          returns in the pointer the new allocation if reallocated.
 * @return the number of bytes of the serialization
 */
size_t serializeHOTChildrenToList_NB(NBodyHOTNode_t* children[8], uint8_t** data_p, size_t* curAlloc);


/**
 * Deserialize a linked list of hashed octree nodes from an array of
 * bytes, data, of size size.
 *
 * @param data: the array of bytes to deserialize
 * @param size: the number of bytes to deserialize
 */
NBodyHOTNode_t* deserializeHOTNodeList_NB(uint8_t* data, size_t size);


/**
 * Serialize an array of 8 hashed octree node children.
 * Some entries in the array may be NULL, but they
 * are still serialized into null-data.
 * The resulting serialization is exactly HOTChildrenSerSize bytes.
 *
 * @param children: the array of HOT node children
 * @param data[out]: an array of bytes of size at least HOTChildrenSerSize
 */
void serializeHOTChildren_NB(NBodyHOTNode_t* children[8], uint8_t* data);


/**
 * Deserialize an array of 8 hashed octree node children.
 * Resulting node pointers are stored in the children array.
 * Some resulting entries in the array may be NULL.
 * The data to deserialize is exactly HOTChildrenSerSize bytes.
 *
 * @param children[out]: the array of HOT node children
 * @param data: the array of bytes encoding the octree children.
 */
void deserializeHOTChildren_NB(NBodyHOTNode_t* children[8], uint8_t* data);


/**
 * Serialize an Octree, tree. Returns the array of bytes
 * in data_p which is the serialization of the tree.
 *
 * @param tree: the tree to serialize
 * @param[in,out] data_p: a pointer to the array of bytes of tree data
 * @param prevAlloc: the current allocation of the array pointed to by data_p
 * @return the size of the data in bytes
 */
size_t serializeOctree_NB(NBOctree_t* tree, uint8_t** data_p, size_t* prevAlloc);


/**
 * Deserialize an Octree from an array of bytes, pointed to by
 * data and of size size.
 *
 * @param data: the array of bytes encoding an Octree
 * @param size: the sizeo the data array
 * @return the deserialized Octree.
 */
NBOctree_t* deserializeOctree_NB(uint8_t* data, size_t size);


#ifdef __cplusplus
}
#endif

#endif
