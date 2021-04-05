
#ifndef _NBODY_KEYS_
#define _NBODY_KEYS_

#include <string.h>
#include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned long long spatialKey_t;

/**
 * The maximum length of one dimenstion in the key in bits.
 */
static const int dimKeyLen = 21;

/**
 * An upper bound on the value of one dimension in the key.
 */
static const spatialKey_t dimKeyMax = (1ull << 21);

/**
 * An upper bound on the value of one dimension in the key as a double.
 */
static const double dimKeyMax_d = (double) (1ull << 21);

/**
 * A mask for extracting the essential part one dimension's intermediate key.
 */
static const spatialKey_t dimKeyMask = 0xfffff80000000000;

/**
 * The shift to apply after applying the mask to normalize a dimension's
 * intermediate key.
 */
static const int dimKeyShift = 43;

/**
 * Masks for the morton bit swizzle.
 */
static const spatialKey_t sepMasks[] = {
	0x1f00000000ffff,
	0x1f0000ff0000ff,
	0x100f00f00f00f00f,
	0x10c30c30c30c30c3,
	0x1249249249249249
};


/**
 * Allocate an array of N spatial keys.
 *
 * @param N: the number of keys to alloc
 * @param[out] keys: a pointer to the created keys.
 * @return 0 iff the allocation was successful
 */
static inline int allocSpatialKeys_NB(long N, spatialKey_t** keys) {
	if (keys == NULL) {
		return 1;
	}

	*keys = (spatialKey_t*) malloc(sizeof(spatialKey_t)*N);
	return (*keys == NULL);
}


static inline int reallocSpatialKeys_NB(long N, spatialKey_t** keys) {
	if (keys == NULL) {
		return 1;
	}
	*keys = (spatialKey_t*) realloc(*keys, sizeof(spatialKey_t)*N);
	return (*keys == NULL);
}

/**
 * Free an array of N spatial keys.
 *
 * @param keys the array to free.
 */
static inline void freeSpatialKeys_NB(spatialKey_t* keys) {
	free(keys);
}


/**
 * Compute the spatial key for position i and store it in keys[i].
 * This is the 3D Morton ordering for X > Y > Z.
 *
 * @param i: the index of the position to compute a key for.
 * @param r: the array positions
 * @param keys: the array of keys to store the computed key in.
 * @param domainSize: the overall size of the space containing all positions.
 */
void computeSpatialKey_NB(long i, const double* r, spatialKey_t* keys, double domainSize);


/**
 * Compute the spatial key for N positions and store then in keys
 * with the corresponding index as the position.
 *
 * @param N: the number of positions and keys.
 * @param r: a 3*N array positions
 * @param keys: the array of keys to store the computed keys in.
 */
void computeSpatialKeys_NB(long N, double* r, spatialKey_t* keys);


/**
 * Compute the spatial key for N positions and store then in keys
 * with the corresponding index as the position.
 * Positions should be bounded by the box [-domainSize, domainSize]
 *
 * @param N: the number of positions and keys.
 * @param domainSize: the ceiling of absolute value of any position in any dimension.
 * @param r: a 3*N array positions
 * @param keys: the array of keys to store the computed keys in.
 */
static inline void computeSpatialKeysDomain_NB(long N, double domainSize, double* r, spatialKey_t* keys) {
	long i;
	for (i = 0; i < N; ++i) {
		computeSpatialKey_NB(i, r, keys, domainSize);
	}
}


/**
 * Sort the N spatial keys as if they are integers.
 * Based on the key ordering also sort the positions and velocities.
 * This function requires pre-allocated space idx, of size N.
 *
 * @param N: the number of keys to sort
 * @param keys: an array of N spatial keys where key[i] is position i.
 * @param idx[out]: an array of size N used as a mapping for parallel arrays
                    whereby v[i] should be moved to i.
 */
void sortSpatialKeys_NB(long N, spatialKey_t* keys, long* idx);


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of 4*N float values.
 * Note that this sorting is done in-place. As a by-product,
 * the idx array contents are also shuffeled in the process.
 * To preserve the index map, a copy should be made before calling
 * this function.
 *
 * @param N: the number of values.
 * @param idx: the index map.
 * @param vals: the 4*N array of float values to sort.
 */
void sortByIdxMap4N_NB(long N, long* idx, float* vals);


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of 3*N double values.
 * Note that this sorting is done in-place. As a by-product,
 * the idx array contents are also shuffeled in the process.
 * To preserve the index map, a copy should be made before calling
 * this function.
 *
 * @param N: the number of values.
 * @param idx: the index map.
 * @param vals: the 3*N array of double values to sort.
 */
void sortByIdxMap3N_NB(long N, long* idx, double* vals);


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of N double values.
 * Note that this sorting is done in-place. As a by-product,
 * the idx array contents are also shuffeled in the process.
 * To preserve the index map, a copy should be made before calling
 * this function.
 *
 * @param N: the number of values.
 * @param idx: the index map.
 * @param vals: the N array of float values to sort.
 */
void sortByIdxMap_NB(long N, long* idx, double* vals);


#if NBODY_MPI
/**
 * Compute spatial keys (morton ordering) of bodies and sort them based on this.
 * The sort is carried out by a distributed sorting routine modified from bubble sort:
 *   - each processor sorts their local data.
 *   - ranks share a fixed number of their highest values with the processor on the right
 *     (and lowest values with rank on the left)
 *   - Highest values on rank i are sorted alonside lowest values from rank i+1
 *   - Rank i keeps lower half of exchanged and sorted elements, rank i+1 keeps higher half.
 *   - Repeat until no changes
 *
 * @note for any rank, localN must be at least twice the SORT_EXCHANGE_SIZE.
 *
 * @param N: total number across all procs, used for optional verification
 * @param localN: number of bodies on local processor
 * @param domainSize: (ceiling of) largest extent in any dimension encompassing all positions
 * @param r, v, m, work: arrays of size 3*localN, 3*localN, localN, localN, holding the data to sort
 * @param keys, idx, idx2: arrays of size localN used as temporary stoarge during the algorithm
 */
void computeAndSortKeys_NBMPI(long N, long localN, double domainSize, spatialKey_t* keys,
    double* r,  double* v,  double* m,  double* work, long* idx, long* idx2);

#endif



#ifdef __cplusplus
}
#endif

#endif
