
#include "NBodyKeys.h"

#include "NBodyHelpers.h"

#include "NBodyConfig.h"

#if NBODY_MPI
#include <mpi.h>
#endif



/**
 * Compute the spatial key for position i and store it in keys[i].
 *
 * @param i, the index of the position to compute a key for.
 * @param r, the array positions
 * @param keys, the array of keys to store the computed key in.
 * @param domainSize, the overall size of the space containing all positions.
 */
void computeSpatialKey_NB(long i, const double* r, spatialKey_t* keys, double domainSize) {

	// map (-R,R) -> (0,2R) -> (0, dimKeyMax)
	// domainSize = R is already ceil of any position in any dimension.
	// Originally, Warren and Salmon of Hashed Octree did
	// (-R,R) -> "integers" -> take (key_len - 1)/numDims most sig. bits
	// Going (0,2R) -> (0, dimKeyMax) is same as taking most sig. bits

	double X[3] = {
		r[3*i]   + domainSize,
		r[3*i+1] + domainSize,
		r[3*i+2] + domainSize
	};

	double frac = dimKeyMax_d / (2.0*domainSize);

	X[0] *= frac;
	X[1] *= frac;
	X[2] *= frac;

	//convert floats to ints.
	spatialKey_t Y[3] = {
		(spatialKey_t) (X[0]),
		(spatialKey_t) (X[1]),
		(spatialKey_t) (X[2])
	};

	//morton code swizzle
	Y[0] = (Y[0] | Y[0] << 32) & sepMasks[0];
	Y[0] = (Y[0] | Y[0] << 16) & sepMasks[1];
	Y[0] = (Y[0] | Y[0] << 8)  & sepMasks[2];
	Y[0] = (Y[0] | Y[0] << 4)  & sepMasks[3];
	Y[0] = (Y[0] | Y[0] << 2)  & sepMasks[4];

	Y[1] = (Y[1] | Y[1] << 32) & sepMasks[0];
	Y[1] = (Y[1] | Y[1] << 16) & sepMasks[1];
	Y[1] = (Y[1] | Y[1] << 8)  & sepMasks[2];
	Y[1] = (Y[1] | Y[1] << 4)  & sepMasks[3];
	Y[1] = (Y[1] | Y[1] << 2)  & sepMasks[4];

	Y[2] = (Y[2] | Y[2] << 32) & sepMasks[0];
	Y[2] = (Y[2] | Y[2] << 16) & sepMasks[1];
	Y[2] = (Y[2] | Y[2] << 8)  & sepMasks[2];
	Y[2] = (Y[2] | Y[2] << 4)  & sepMasks[3];
	Y[2] = (Y[2] | Y[2] << 2)  & sepMasks[4];

	//prepend 1, bit-wise or to interleave.
	keys[i] = (1ull << 63) | (Y[0] << 2) | (Y[1] << 1) | Y[2];
}


void computeSpatialKeys_NB(long N, double* r, spatialKey_t* keys) {
	double domainSize = computeDomainSize_NB(N, r);
	long i;
	for (i = 0; i < N; ++i) {
		computeSpatialKey_NB(i, r, keys, domainSize);
	}
}


/**
 * Sort the N spatial keys as if they are integers.
 * Based on the key ordering also sort the positions and velocities.
 * This function requires pre-allocated space idx and idx2, each of size N.
 *
 * @param N, the number of keys to sort.
 * @param keys, an array of N spatial keys where key[i] is position i.
 * @param idx[out], a mapping of how the sorting occured. idx[i] should move to i.
 */
void sortSpatialKeys_NB(long N, spatialKey_t* keys, long* idx) {
//note: don't need to sort accels since they are overwritten completely in next step

	long i;
	//Indices store where values should be eventually moved by sorting of keys.
	//idx[i] should move to i.
	//Do this instead of sorting keys, r, and v simultaneously for better locality.
	for (i = 0; i < N; ++i) {
		idx[i] = i;
	}

	//sort by insertion sort since vals should be in approximately
	//correct order by prev step's sorting.
	long j, curIdx;
	spatialKey_t ki;
	for (i = 0; i < N; ++i) {
		ki = keys[i];
		curIdx = idx[i];
		for (j = i-1; j >= 0 && keys[j] > ki; --j) {
			keys[j+1] = keys[j];
			idx[j+1] = idx[j];
		}
		idx[j+1] = curIdx;
		keys[j+1] = ki;
	}
}


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of 3*N double values.
 *
 * @param N, the number of values.
 * @param idx, the index map.
 * @param vals, the 3*N array of double values to sort.
 */
void sortByIdxMap3N_NB(long N, long* idx, double* vals) {


	//sort r's; recall idx[i] should move to i.
	double tmp[3];
	long src, dst;
	for (long i = 0; i < N; ++i) {
		dst = i;
		while (idx[dst] != i) {
			src = idx[dst];

			//swap src and dst;
			tmp[0] = vals[3*src];
			tmp[1] = vals[3*src + 1];
			tmp[2] = vals[3*src + 2];

			vals[3*src]     = vals[3*dst];
			vals[3*src + 1] = vals[3*dst + 1];
			vals[3*src + 2] = vals[3*dst + 2];

			vals[3*dst]     = tmp[0];
			vals[3*dst + 1] = tmp[1];
			vals[3*dst + 2] = tmp[2];

			//update idx to say dst is properly placed
			idx[dst] = dst;
			//continue the train of swaps.
			//New dst is whereever was just swapped to, new src is based on idx.
			dst = src;
			src = idx[dst];
		}

		//if we broke, then the swapping was successful and record so.
		idx[dst] = dst;
	}

}


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of 4*N float values.
 *
 * @param N, the number of values.
 * @param idx, the index map.
 * @param vals, the 4*N array of float values to sort.
 */
void sortByIdxMap4N_NB(long N, long* idx, float* vals) {

	//sort r's; recall idx[i] should move to i.
	float tmp[4];
	long src, dst;
	for (long i = 0; i < N; ++i) {
		dst = i;
		while (idx[dst] != i) {
			src = idx[dst];

			//swap src and dst;
			tmp[0] = vals[4*src];
			tmp[1] = vals[4*src + 1];
			tmp[2] = vals[4*src + 2];
			tmp[3] = vals[4*src + 3];

			vals[4*src]     = vals[4*dst];
			vals[4*src + 1] = vals[4*dst + 1];
			vals[4*src + 2] = vals[4*dst + 2];
			vals[4*src + 3] = vals[4*dst + 3];

			vals[4*dst]     = tmp[0];
			vals[4*dst + 1] = tmp[1];
			vals[4*dst + 2] = tmp[2];
			vals[4*dst + 3] = tmp[3];

			//update idx to say dst is properly placed
			idx[dst] = dst;
			//continue the train of swaps.
			//New dst is whereever was just swapped to, new src is based on idx.
			dst = src;
			src = idx[dst];
		}

		//if we broke, then the swapping was successful and record so.
		idx[dst] = dst;
	}

}


/**
 * Based on an index map where index idx[i] should move to index i,
 * sort an array of N double values.
 *
 * @param N, the number of values.
 * @param idx, the index map.
 * @param vals, the N array of float values to sort.
 */
void sortByIdxMap_NB(long N, long* idx, double* vals) {

	//sort r's; recall idx[i] should move to i.
	double tmp[3];
	long src, dst;
	for (long i = 0; i < N; ++i) {
		dst = i;
		while (idx[dst] != i) {
			src = idx[dst];

			tmp[0] = vals[src];
			vals[src] = vals[dst];
			vals[dst] = tmp[0];

			idx[dst] = dst;
			dst = src;
			src = idx[dst];
		}

		idx[dst] = dst;
	}

}


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
    double* r,  double* v,  double* m,  double* work, long* idx, long* idx2)
{

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    computeSpatialKeysDomain_NB(localN, domainSize, r, keys);

    spatialKey_t* exchangeKeysLeft = (spatialKey_t*) malloc(sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE*2);
    spatialKey_t* exchangeKeysRight = (spatialKey_t*) malloc(sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE*2);
    double* exchangeRLeft = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2*3);
    double* exchangeRRight = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2*3);
    double* exchangeVLeft = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2*3);
    double* exchangeVRight = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2*3);
    double* exchangeMLeft = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2);
    double* exchangeMRight = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2);
    double* exchangeWorkLeft = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2);
    double* exchangeWorkRight = (double*) malloc(sizeof(double)*SORT_EXCHANGE_SIZE*2);

    /*
     * Here the idea is for adjacent processors (assume line topology of increasing rank)
     * to exchange a set of elements, both sort the exchanged elements,
     * and then keep their respective halves (lower rank takes lower half, and vice versa).
     * E.g., Rank 0 sends highest SORT_EXCHANGE_SIZE elements to rank 1, 1 sends lowest SORT_EXCHANGE_SIZE
     * elements to 0. This 2*SORT_EXCHANGE_SIZE array is sorted. 0 keeps the lower half,
     * discarding the top half. 1 keeps the higher half, discarding the bottom half.
     * All non-boundary processors (i.e. not 0 or world_size-1) have two adjacent ranks.,
     * Even processors send right first, then left second. Odd processors are opposite.
     * Repeat until there is no change.
     */

    // int steps = 0;
    while(1) {
        // steps++;
        sortSpatialKeys_NB(localN, keys, idx);
        memcpy(idx2, idx, sizeof(long)*localN);
        sortByIdxMap3N_NB(localN, idx2, r);
        memcpy(idx2, idx, sizeof(long)*localN);
        sortByIdxMap3N_NB(localN, idx2, v);
        memcpy(idx2, idx, sizeof(long)*localN);
        sortByIdxMap_NB(localN, idx2, m);
        sortByIdxMap_NB(localN, idx, work);


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

        //copy data into exchange buffers
        memcpy(exchangeKeysLeft, keys, sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeKeysLeft + SORT_EXCHANGE_SIZE, exchangeKeysLeft, sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeRLeft, r, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeRLeft + 3*SORT_EXCHANGE_SIZE, exchangeRLeft, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeVLeft, v, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeVLeft + 3*SORT_EXCHANGE_SIZE, exchangeVLeft, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeMLeft, m, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeMLeft + SORT_EXCHANGE_SIZE, exchangeMLeft, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeWorkLeft, work, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeWorkLeft + SORT_EXCHANGE_SIZE, exchangeWorkLeft, sizeof(double)*SORT_EXCHANGE_SIZE);

        memcpy(exchangeKeysRight, keys + (localN - SORT_EXCHANGE_SIZE), sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeKeysRight + SORT_EXCHANGE_SIZE, exchangeKeysRight, sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeRRight, r + 3*(localN - SORT_EXCHANGE_SIZE), sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeRRight + 3*SORT_EXCHANGE_SIZE, exchangeRRight, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeVRight, v + 3*(localN - SORT_EXCHANGE_SIZE), sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeVRight + 3*SORT_EXCHANGE_SIZE, exchangeVRight, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(exchangeMRight, m + (localN - SORT_EXCHANGE_SIZE), sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeMRight + SORT_EXCHANGE_SIZE, exchangeMRight, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeWorkRight, work + (localN - SORT_EXCHANGE_SIZE), sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(exchangeWorkRight + SORT_EXCHANGE_SIZE, exchangeWorkRight, sizeof(double)*SORT_EXCHANGE_SIZE);


        //exchange right and then left
        MPI_Status stat;
        if (firstExchangePartner >= 0) {
            if (firstBuffRight) {
                MPI_Sendrecv_replace(exchangeKeysRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_UNSIGNED_LONG_LONG,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeRRight + 3*SORT_EXCHANGE_SIZE, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeVRight + 3*SORT_EXCHANGE_SIZE, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeMRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeWorkRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
            } else {
                MPI_Sendrecv_replace(exchangeKeysLeft, SORT_EXCHANGE_SIZE, MPI_UNSIGNED_LONG_LONG,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeRLeft, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeVLeft, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeMLeft, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeWorkLeft, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     firstExchangePartner, 0, firstExchangePartner, 0, MPI_COMM_WORLD, &stat);
            }
        }
        if (secondExchangePartner >= 0) {
            if (firstBuffRight) {
                MPI_Sendrecv_replace(exchangeKeysLeft, SORT_EXCHANGE_SIZE, MPI_UNSIGNED_LONG_LONG,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeRLeft, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeVLeft, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeMLeft, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeWorkLeft, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
            } else {
                MPI_Sendrecv_replace(exchangeKeysRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_UNSIGNED_LONG_LONG,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeRRight + 3*SORT_EXCHANGE_SIZE, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeVRight + 3*SORT_EXCHANGE_SIZE, 3*SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeMRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv_replace(exchangeWorkRight + SORT_EXCHANGE_SIZE, SORT_EXCHANGE_SIZE, MPI_DOUBLE,
                                     secondExchangePartner, 0, secondExchangePartner, 0, MPI_COMM_WORLD, &stat);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //if there was an exchange, sort the keys and then sort remaining arrays by index.
        if (firstExchangePartner >= 0) {
            if (world_rank % 2 == 0) {
                sortSpatialKeys_NB(2*SORT_EXCHANGE_SIZE, exchangeKeysRight, idx);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeRRight);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeVRight);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeMRight);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx, exchangeWorkRight);

            } else {
                sortSpatialKeys_NB(2*SORT_EXCHANGE_SIZE, exchangeKeysLeft, idx);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeRLeft);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeVLeft);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeMLeft);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx, exchangeWorkLeft);
            }
        }
        if (secondExchangePartner >= 0) {
            if (world_rank % 2 == 0) {
                sortSpatialKeys_NB(2*SORT_EXCHANGE_SIZE, exchangeKeysLeft, idx);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeRLeft);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeVLeft);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeMLeft);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx, exchangeWorkLeft);
            } else {
                sortSpatialKeys_NB(2*SORT_EXCHANGE_SIZE, exchangeKeysRight, idx);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeRRight);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap3N_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeVRight);
                memcpy(idx2, idx, sizeof(long)*2*SORT_EXCHANGE_SIZE);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx2, exchangeMRight);
                sortByIdxMap_NB(2*SORT_EXCHANGE_SIZE, idx, exchangeWorkRight);
            }
        }

        int changed = 0;
        //compare against top half of leftExchange for changes
        for (int i = 0; i < SORT_EXCHANGE_SIZE && !changed; ++i) {
            if (keys[i] != exchangeKeysLeft[SORT_EXCHANGE_SIZE + i]) {
                changed = 1;
            }
        }
        //compare against bottom half of rightExchange for changes
        for (int i = 1; i <= SORT_EXCHANGE_SIZE && !changed; ++i) {
            if (keys[localN - i] != exchangeKeysRight[SORT_EXCHANGE_SIZE - i]) {
                changed = 1;
            }
        }

        //if no changes for all procs, then we are done sorting
        int anyChanged = 0;
        MPI_Allreduce(&changed, &anyChanged, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
        if (anyChanged == 0) {
            break;
        }

        //keep top-half of left neighbour group, bottom-half of right neighbour group
        memcpy(keys, exchangeKeysLeft + SORT_EXCHANGE_SIZE, sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(keys + (localN - SORT_EXCHANGE_SIZE), exchangeKeysRight, sizeof(spatialKey_t)*SORT_EXCHANGE_SIZE);
        memcpy(r, exchangeRLeft + 3*SORT_EXCHANGE_SIZE, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(r + 3*(localN - SORT_EXCHANGE_SIZE), exchangeRRight, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(v, exchangeVLeft + 3*SORT_EXCHANGE_SIZE, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(v + 3*(localN - SORT_EXCHANGE_SIZE), exchangeVRight, sizeof(double)*SORT_EXCHANGE_SIZE*3);
        memcpy(m, exchangeMLeft + SORT_EXCHANGE_SIZE, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(m + (localN - SORT_EXCHANGE_SIZE), exchangeMRight, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(work, exchangeWorkLeft + SORT_EXCHANGE_SIZE, sizeof(double)*SORT_EXCHANGE_SIZE);
        memcpy(work + (localN - SORT_EXCHANGE_SIZE), exchangeWorkRight, sizeof(double)*SORT_EXCHANGE_SIZE);
    }

    free(exchangeKeysLeft);
    free(exchangeKeysRight);
    free(exchangeRRight);
    free(exchangeRLeft);
    free(exchangeVRight);
    free(exchangeVLeft);
    free(exchangeMLeft);
    free(exchangeMRight);
    free(exchangeWorkRight);
    free(exchangeWorkLeft);

#if 0
    //verify distributed sort by gather and local sort
    spatialKey_t* allKeys = NULL;
    if (world_rank == 0) {
        allKeys = (spatialKey_t*) malloc(sizeof(spatialKey_t)*N);
    }
    MPI_Gather(keys, localN, MPI_UNSIGNED_LONG_LONG,
                allKeys, localN, MPI_UNSIGNED_LONG_LONG,
                0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        fprintf(stderr, "allKeys = [");
        for (int i = 0; i < N; ++i) {
            fprintf(stderr, "%llu ", allKeys[i] );
        }
        fprintf(stderr, "]\n\n");

        spatialKey_t* sortedKeys = (spatialKey_t*) malloc(sizeof(spatialKey_t)*N);
        memcpy(sortedKeys, allKeys, sizeof(spatialKey_t)*N);

        long* tmpIdx = (long*) malloc(sizeof(long)*N);
        sortSpatialKeys_NB(N, sortedKeys, tmpIdx);

        fprintf(stderr, "sortedKeys = [");
        for (int i = 0; i < N; ++i) {
            fprintf(stderr, "%llu ", sortedKeys[i] );
        }
        fprintf(stderr, "]\n\n");

        for (int i = 0; i < N; ++i) {
            if (sortedKeys[i] != allKeys[i]) {
                fprintf(stderr, "LISTS WERE NOT EQUAL AT IDX %d: %llu vs %llu\n", i, sortedKeys[i], allKeys[i]);
            }
        }

        free(sortedKeys);
        free(tmpIdx);
        free(allKeys);
        fprintf(stderr, "\nDONE VERIFY.\n");

    }
#endif
}
#endif