
#include "NBodyHelpers.h"

#if NBODY_MPI
#include <mpi.h>

/**
 * Smartly distribute the bodies (i.e. load balance)
 * across all processors so that work is shared evenly.
 * Work estimates from previous simulation step used. Assumes bodies and additional data
 * are already sorted and now the processor boundaries must be adjusted to balance the work.
 *
 * This method does not do a gather/scatter. Rather, data is scanned left to right
 * from rank 0 to rank world_size-1, shifting data on processor boundaries
 * left and right to enforce a load balance.
 *
 * Data is passed as double pointers in case of reallocation.
 *
 * @param N: total number of bodies across all processors
 * @param localN_p: a pointer to the number of bodies local to this processor
 * @param currAlloc_p: a pointer to the current allocation size of each data array,
                       (or 1/3 the actual allocation size if the array is a triple, i.e. pos and velocity)
 * @param r, v, a: arrays of pos, velocity, accel data of size 3*localN
 * @param m: array of mass of size localN
 * @param work: array of work estimates of size localN
 * @param keys, idx1, idx2: temporary space of allocation *currAlloc_p whose allocation may change
 */
void distributeWork_NBMPI(long N, long* localN_p, long* currAlloc_p,
    double** r, double** v, double** a, double** m, double** work,
    spatialKey_t** keys, long** idx1, long** idx2) {

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    long localN = *localN_p;
    long currAlloc = *currAlloc_p;

    //compute total work and per-proc work target
    double* w = *work;
    double totalWork = 0.0;
    for (int i = 0; i < localN; ++i) {
        totalWork += w[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &totalWork, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double workTarg = totalWork / world_size;

    //iterate to compute current "cost zone"
    //rank world_size-1 just gets left with everything left over
    MPI_Status stat;
    for (int k = 0; k < world_size-1; ++k) {
        if (world_rank == k) {
            totalWork = 0.0;
            long sendIdx = -1;
            w = *work;
            for (long i = 0; i < localN; ++i) {
                totalWork += w[i];
                if (totalWork >= workTarg && sendIdx == -1) {
                    //then data from 0..i, inclusive should be kept by this proc.
                    //don't break so we compute total work for workDiff
                    sendIdx = i + 1;
                }
            }

            //compute workDiff which represents how much additional work this proc has
            //w.r.t. the work target. a negative value means we want to take some data
            //from the right-neighbour to make up the difference. A positive number
            //means we have too much and want to push data to the right neighbour
            double workDiff = totalWork - workTarg;
            if (sendIdx == localN) {
                workDiff = 0.0;
            }
            MPI_Send(&workDiff, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);

            if (workDiff < 0.0) {
                //pull data from right neighbour
                int recvCount = 0;
                MPI_Recv(&recvCount, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, &stat);
		if (localN + recvCount > currAlloc) {
                    while (localN + recvCount > currAlloc) {
                        currAlloc += (N / world_size);
                    }
                    int err;
                    err = reallocData_NB(currAlloc, r, v, a, m, work);
                    err |= reallocSpatialKeys_NB(currAlloc, keys);
                    err |= !(*idx1 = (long*) realloc(*idx1, sizeof(long)*currAlloc));
                    err |= !(*idx2 = (long*) realloc(*idx2, sizeof(long)*currAlloc));
                    if (err) {
                        fprintf(stderr, "Could not alloc local data for N=%ld\n", currAlloc);
                        exit(ALLOC_ERROR);
                    }
                }
                MPI_Recv((*r) + 3*localN, 3*recvCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv((*v) + 3*localN, 3*recvCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv((*m) + localN, recvCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv((*work) + localN, recvCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, &stat);
                localN = localN + recvCount;
            } else if (workDiff > 0.0) {
                //push data to right neighbour
                int sendCount = localN - sendIdx;
                MPI_Send(&sendCount, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send((*r) + 3*sendIdx, 3*sendCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send((*v) + 3*sendIdx, 3*sendCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send((*m) + sendIdx, sendCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send((*work) + sendIdx, sendCount, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
                localN = localN - sendCount;
            }
        } else if (world_rank == k + 1) {
            double workDiff = 0.0;
            MPI_Recv(&workDiff, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
	    if (workDiff < 0.0) {
                //push data to prev processor
                workDiff *= -1.0;
                totalWork = 0.0;
                w = *work;
                int sendCount = -1;
                for (int i = 0; i < localN; ++i) {
                    totalWork += w[i];
                    if (totalWork >= workDiff) {
                        sendCount = i+1;
                        break;
                    }
                }
                if (sendCount < 0) {
                    sendCount = localN;
                }
                MPI_Send(&sendCount, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD);
                MPI_Send(*r, 3*sendCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
                MPI_Send(*v, 3*sendCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
                MPI_Send(*m, sendCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
                MPI_Send(*work, sendCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);

                int newN = localN - sendCount;
                memmove(*r, (*r) + 3*sendCount, sizeof(double)*newN*3);
                memmove(*v, (*v) + 3*sendCount, sizeof(double)*newN*3);
                memmove(*m, (*m) + sendCount, sizeof(double)*newN);
                memmove(*work, (*work) + sendCount, sizeof(double)*newN);
                localN = newN;
            } else if (workDiff > 0.0) {
                //pull excess data from prev processor
                int recvCount = 0;
                MPI_Recv(&recvCount, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
                if (localN + recvCount > currAlloc) {
                    while (localN + recvCount > currAlloc) {
                        currAlloc += (N / world_size);
                    }
                    int err;
                    err = reallocData_NB(currAlloc, r, v, a, m, work);
                    err |= reallocSpatialKeys_NB(currAlloc, keys);
                    err |= !(*idx1 = (long*) realloc(*idx1, sizeof(long)*currAlloc));
                    err |= !(*idx2 = (long*) realloc(*idx2, sizeof(long)*currAlloc));
                    if (err) {
                        fprintf(stderr, "Could not alloc local data for N=%ld\n", currAlloc);
                        exit(ALLOC_ERROR);
                    }
                }

                //shift current data right to make room for extra;
                memmove((*r) + 3*recvCount, *r, sizeof(double)*localN*3);
                memmove((*v) + 3*recvCount, *v, sizeof(double)*localN*3);
                memmove((*m) + recvCount, *m, sizeof(double)*localN);
                memmove((*work) + recvCount, *work, sizeof(double)*localN);

                MPI_Recv(*r, 3*recvCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv(*v, 3*recvCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv(*m, recvCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv(*work, recvCount, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, &stat);
                
                localN = localN + recvCount;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    *localN_p = localN;
    *currAlloc_p = currAlloc;
}

#endif


