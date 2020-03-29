
#include "NBodyKeys.h"

/**
 * Compute the spatial key for position i and store it in keys[i].
 *
 * @param i, the index of the position to compute a key for.
 * @param r, the array positions
 * @param keys, the array of keys to store the computed key in.
 * @param domainSize, the overall size of the space containing all positions.
 */
void computeSpatialKey_NB(long i, double* r, spatialKey_t* keys, double domainSize) {

	// map [-R,R] -> [0,2R] -> [0, dimKeyMax)
	// Originally, authors of Hashed Octree did
	// (-R,R) -> "integers" -> take (key_len - 1)/numDims most sig. bits

	double X[3] = {
		r[3*i]   + domainSize,
		r[3*i+1] + domainSize,
		r[3*i+2] + domainSize
	};

	//A little extra to ensure X*frac in [0,dimKeyMax)
	double frac = dimKeyMax_d / (2.00001 * domainSize);

	X[0] *= frac;
	X[1] *= frac;
	X[2] *= frac;

	// fprintf(stderr, "X[0]: %g\n", ceil(X[0]));
	// fprintf(stderr, "X[1]: %g\n", ceil(X[1]));
	// fprintf(stderr, "X[2]: %g\n", ceil(X[2]));

	spatialKey_t Y[3] = {
		(spatialKey_t) ceil(X[0]),
		(spatialKey_t) ceil(X[1]),
		(spatialKey_t) ceil(X[2])
	};

	// fprintf(stderr, "Y[0]: %llx\n", Y[0]);
	// fprintf(stderr, "Y[1]: %llx\n", Y[1]);
	// fprintf(stderr, "Y[2]: %llx\n", Y[2]);

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
	long src, dst, i;
	for (i = 0; i < N; ++i) {
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
	long src, dst, i;
	for (i = 0; i < N; ++i) {
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
	long src, dst, i;
	for (i = 0; i < N; ++i) {
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
