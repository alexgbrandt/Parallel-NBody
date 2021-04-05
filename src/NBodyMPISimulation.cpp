
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "NBodyMPISimulation.hpp"

#include "NBodyConfig.h"

#include "NBodyInit.h"
#include "NBodyHelpers.h"
#include "NBodyForces.h"
#include "NBodyIntegrator.h"

#include "NBodyKeys.h"
#include "NBodyOctree.h"
#include "NBodyParallel.h"

#if NBODY_SIM_WITH_RENDERER
#include "ogl/NBodyRenderer.hpp"
#endif

#include <unistd.h>
#include <mpi.h>
#include "NBodySerialize.h"

#include "Unix_Timer.h"


void NBodyNaive_1(long N, double dt, double t_end, time_t seed) {

/**********************************
 * Init data
 **********************************/
    int err;
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL;
    err = allocData3N_NB(N, &r);
    err |= allocData3N_NB(N, &v);
    err |= allocData3N_NB(N, &a);
    err |= allocDataN_NB(N, &m);
    if (err) {
        fprintf(stderr, "Could not alloc data for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }

    err = initData_NB(N, seed, r, v, a, m);

    double Epot = computeEpot_NB(N, m, r);
    double Ekin = computeEkin_NB(N, m, v);
    double E0 = Epot + Ekin;
    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "E0: %.15g\n", E0);


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = createColors(N);
    NBodyRenderer renderer;
    renderer.startRenderThread();
    // renderer.updatePoints(N, r);

    // renderer.setOctree(octree);
    // renderer.joinRenderThread();
    // return 0;
#endif


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;
    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
#if NBODY_SIM_WITH_RENDERER
        if (renderer.shouldClose()) {
            break;
        }
#endif
        //kick and drift
        performNBodyHalfStepA(N, dt, r, v, a, m);

#if NBODY_SIM_WITH_RENDERER
        //update renderer before next kick since positions only change in step A
        if (renderer.needsUpdate()) {
            renderer.updatePoints(N, r, m, 1, colors);
        }
#endif

        //update accels
        computeForces(N, m, r, a);

        //second kick
        performNBodyHalfStepB(N, dt, r, v, a, m);

        // Print checkpoints.
        // if (t >= t_out) {
        //     Ekin = computeEkin_NB(N, m, v);
        //     Epot = computeEpot_NB(N, m, r);
        //     double E1 = Ekin + Epot;
        //     for (long k = 0; k < N; ++k) {
        //         printf("%.15g %.15g %.15g ", r[3*k], r[3*k+1], r[3*k+2]);
        //     }
        //     printf("%.15g\n", E1);

        //     t_out += dt_out;
        // }
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    Epot = computeEpot_NB(N, m, r);
    Ekin = computeEkin_NB(N, m, v);
    E0 = Epot + Ekin;
    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "Eend: %.15g\n", E0);
    fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);


#if NBODY_SIM_WITH_RENDERER
    //End of simulation, wait for renderer to finish.
    renderer.joinRenderThread();
    //End of simultation, force renderer to finish.
    // renderer.stopRenderThread();
#endif


/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
}



void NBodyNaive_NBMPI(long N, double dt, double t_end, time_t seed) {

/**********************************
 * Init data
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL;
    double *r_rec = NULL, *m_rec = NULL;

    //for local proc's bodies
    err = allocData3N_NB(localN, &r);
    err = allocData3N_NB(localN, &v);
    err = allocData3N_NB(localN, &a);
    err = allocDataN_NB(localN, &m);

    //working space for receiving other's data
    err = allocData3N_NB(localN, &r_rec);
    err = allocDataN_NB(localN, &m_rec);

    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData3N_NB(N, &rr);
        err |= allocData3N_NB(N, &vv);
        err |= allocData3N_NB(N, &aa);
        err |= allocDataN_NB(N, &mm);

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }

    err = MPI_Scatter(rr, localN*3, MPI_DOUBLE,
                       r, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(vv, localN*3, MPI_DOUBLE,
                       v, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(aa, localN*3, MPI_DOUBLE,
                       a, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(mm, localN, MPI_DOUBLE,
                       m, localN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;
    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

#if NBODY_SIM_WITH_RENDERER
        MPI_Gather( r, localN*3, MPI_DOUBLE,
                   rr, localN*3, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            //Note that masses doesn't change here because we aren't doing any sorting.
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesNaive_NBMPI(localN, m, r, a, m_rec, r_rec);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Gather( r, localN*3, MPI_DOUBLE,
               rr, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather( v, localN*3, MPI_DOUBLE,
               vv, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(rr);
    free(vv);
    free(aa);
    free(mm);

    MPI_Finalize();
}


//no spatial decomp or load balancing
//gather, build tree, broadcast
void NBodyBH1_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init data
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL;

    //for local proc's bodies
    err = allocData3N_NB(localN, &r);
    err = allocData3N_NB(localN, &v);
    err = allocData3N_NB(localN, &a);
    err = allocDataN_NB(localN, &m);

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData3N_NB(N, &rr);
        err |= allocData3N_NB(N, &vv);
        err |= allocData3N_NB(N, &aa);
        err |= allocDataN_NB(N, &mm);

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }

    err = MPI_Scatter(rr, localN*3, MPI_DOUBLE,
                       r, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(vv, localN*3, MPI_DOUBLE,
                       v, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(aa, localN*3, MPI_DOUBLE,
                       a, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(mm, localN, MPI_DOUBLE,
                       m, localN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        //Gather all positions in order to build octree
        MPI_Gather( r, localN*3, MPI_DOUBLE,
                   rr, localN*3, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            //build Octree and compute forces to update acceleration.
            domainSize = computeDomainSize_NB(N, rr);
            updateInplace = buildOctreeInPlace_NB(N, rr, mm, domainSize, &octree);
            // fprintf(stderr, "FROM ROOT OCTREE IS: \n");
            // printOctree_NB(octree);
            treeSize = serializeOctree_NB(octree, &treeData, &treeAlloc);
            // for (int j = 0; j < treeSize/128; ++j) {
            //     for (int i = 0; i < 16; ++i) {
            //         // fprintf(stderr, BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(treeData[j*16 + i]));
            //         fprintf(stderr, "%lx", *(((long*) treeData) + (16*j) + i));
            //         fprintf(stderr, " ");
            //     }
            //     fprintf(stderr, "\n");
            // }
        }

        MPI_Bcast(&treeSize, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (treeSize > treeAlloc) {
                treeData = (uint8_t*) realloc(treeData, sizeof(uint8_t)*treeSize);
                treeAlloc = treeSize;
            }
        }
        MPI_Bcast(treeData, treeSize, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            // fprintf(stderr, "\nFROM RANK %d OCTREE DATA IS: \n", world_rank);
            // for (int j = 0; j < treeSize/128; ++j) {
            //     for (int i = 0; i < 16; ++i) {
            //         // fprintf(stderr, BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(treeData[j*16 + i]));
            //         fprintf(stderr, "%lx", *(((long*) treeData) + (16*j) + i));
            //         fprintf(stderr, " ");
            //     }
            //     fprintf(stderr, "\n");
            // }
            freeOctree_NB(octree);
            octree = deserializeOctree_NB(treeData, treeSize);
            // fprintf(stderr, "\nFROM RANK %d OCTREE IS: \n", world_rank);
            // printOctree_NB(octree);
        }

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            //Note that masses doesn't change here because we aren't doing any sorting.
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesOctreeBH_NB(localN, m, r, NULL, a, octree, tmp1_node, tmp2_node, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);

    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Gather( r, localN*3, MPI_DOUBLE,
               rr, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather( v, localN*3, MPI_DOUBLE,
               vv, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(tmp1_node);
    free(tmp2_node);

    MPI_Finalize();
}


//no spatial decomp or load balancing
//build local tree, reduce, broadcast
void NBodyBH2_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init data
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL;

    //for local proc's bodies
    err = allocData3N_NB(localN, &r);
    err = allocData3N_NB(localN, &v);
    err = allocData3N_NB(localN, &a);
    err = allocDataN_NB(localN, &m);

    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData3N_NB(N, &rr);
        err |= allocData3N_NB(N, &vv);
        err |= allocData3N_NB(N, &aa);
        err |= allocDataN_NB(N, &mm);

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }

    err = MPI_Scatter(rr, localN*3, MPI_DOUBLE,
                       r, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(vv, localN*3, MPI_DOUBLE,
                       v, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(aa, localN*3, MPI_DOUBLE,
                       a, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(mm, localN, MPI_DOUBLE,
                       m, localN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;



    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);



    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        MPI_Gather( r, localN*3, MPI_DOUBLE,
                   rr, localN*3, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
        if (world_rank == 0) {
            domainSize = computeDomainSize_NB(N, rr);
        }
        MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        mapReduceBuildOctrees_NBMPI(r, m, localN, domainSize, &octree, &treeData, &treeAlloc);

        if (world_rank == 0) {
            treeSize = serializeOctree_NB(octree, &treeData, &treeAlloc);
            // fprintf(stderr, "\nFROM RANK %d OCTREE IS: \n", world_rank);
            // printOctree_NB(octree);
            // fprintf(stderr, "tree size is: %ld\n", treeSize);
        }

        MPI_Bcast(&treeSize, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (treeSize > treeAlloc) {
                treeData = (uint8_t*) realloc(treeData, sizeof(uint8_t)*treeSize);
                treeAlloc = treeSize;
            }
        }
        MPI_Bcast(treeData, treeSize, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            // fprintf(stderr, "\nFROM RANK %d OCTREE DATA IS: \n", world_rank);
            // for (int j = 0; j < treeSize/128; ++j) {
            //     for (int i = 0; i < 16; ++i) {
            //         // fprintf(stderr, BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(treeData[j*16 + i]));
            //         fprintf(stderr, "%lx", *(((long*) treeData) + (16*j) + i));
            //         fprintf(stderr, " ");
            //     }
            //     fprintf(stderr, "\n");
            // }
            freeOctree_NB(octree);
            octree = deserializeOctree_NB(treeData, treeSize);
            // fprintf(stderr, "\nFROM RANK %d OCTREE IS: \n", world_rank);
            // printOctree_NB(octree);
            // fprintf(stderr, "tree size is: %ld\n", treeSize);
        }

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            //Note that masses doesn't change here because we aren't doing any sorting.
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesOctreeBH_NB(localN, m, r, NULL, a, octree, tmp1_node, tmp2_node, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Gather( r, localN*3, MPI_DOUBLE,
               rr, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather( v, localN*3, MPI_DOUBLE,
               vv, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(tmp1_node);
    free(tmp2_node);

    MPI_Finalize();
}


//spatial decomp but no load balancing
//build local tree, reduce, broadcast
void NBodyBH3_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Alloc data
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL;

    //for local proc's bodies
    err = allocData3N_NB(localN, &r);
    err = allocData3N_NB(localN, &v);
    err = allocData3N_NB(localN, &a);
    err = allocDataN_NB(localN, &m);



/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * Init data
 **********************************/

    spatialKey_t* keys = NULL;
    long* idx = NULL, *idx2 = NULL;
    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData3N_NB(N, &rr);
        err |= allocData3N_NB(N, &vv);
        err |= allocData3N_NB(N, &aa);
        err |= allocDataN_NB(N, &mm);
        err |= allocSpatialKeys_NB(N, &keys);

        //if malloc returns null, set err bit
        err |= !(idx = (long*) malloc(sizeof(long)*N));
        err |= !(idx2 = (long*) malloc(sizeof(long)*N));

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);

        //sort before scattering
        computeSpatialKeys_NB(N, rr, keys);
        sortSpatialKeys_NB(N, keys, idx);
#if NBODY_SIM_WITH_RENDERER
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap4N_NB(N, idx2, colors);
#endif
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, idx2, rr);
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, idx2, vv);
        sortByIdxMap_NB(N, idx, mm);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }

    //scatter initialized data
    err = MPI_Scatter(rr, localN*3, MPI_DOUBLE,
                       r, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(vv, localN*3, MPI_DOUBLE,
                       v, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(aa, localN*3, MPI_DOUBLE,
                       a, localN*3, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Scatter(mm, localN, MPI_DOUBLE,
                       m, localN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        //gather for sorting. We don't need accels so don't gather that.
        err = MPI_Gather( r, localN*3, MPI_DOUBLE,
                         rr, localN*3, MPI_DOUBLE,
                          0, MPI_COMM_WORLD);
        err = MPI_Gather( v, localN*3, MPI_DOUBLE,
                         vv, localN*3, MPI_DOUBLE,
                          0, MPI_COMM_WORLD);
        err = MPI_Gather( m, localN, MPI_DOUBLE,
                         mm, localN, MPI_DOUBLE,
                          0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            //compute keys, sort them, then sort each data array by index map
            //this avoids simultaneous sorting SoA's or having a AoS.
            computeSpatialKeys_NB(N, rr, keys);
            sortSpatialKeys_NB(N, keys, idx);
#if NBODY_SIM_WITH_RENDERER
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap4N_NB(N, idx2, colors);
#endif
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, rr);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, vv);
            sortByIdxMap_NB(N, idx, mm);
        }

        err = MPI_Scatter(rr, localN*3, MPI_DOUBLE,
                           r, localN*3, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
        err = MPI_Scatter(vv, localN*3, MPI_DOUBLE,
                           v, localN*3, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
        err = MPI_Scatter(mm, localN, MPI_DOUBLE,
                           m, localN, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            domainSize = computeDomainSize_NB(N, rr);
        }
        MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        mapReduceBuildOctrees_NBMPI(r, m, localN, domainSize, &octree, &treeData, &treeAlloc);

        if (world_rank == 0) {
            treeSize = serializeOctree_NB(octree, &treeData, &treeAlloc);
        }

        MPI_Bcast(&treeSize, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (treeSize > treeAlloc) {
                treeData = (uint8_t*) realloc(treeData, sizeof(uint8_t)*treeSize);
                treeAlloc = treeSize;
            }
        }
        MPI_Bcast(treeData, treeSize, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            freeOctree_NB(octree);
            octree = deserializeOctree_NB(treeData, treeSize);
        }

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesOctreeBH_NB(localN, m, r, NULL, a, octree, tmp1_node, tmp2_node, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Gather( r, localN*3, MPI_DOUBLE,
               rr, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather( v, localN*3, MPI_DOUBLE,
               vv, localN*3, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gather( m, localN, MPI_DOUBLE,
               mm, localN, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(tmp1_node);
    free(tmp2_node);
    free(idx);
    free(idx2);
    free(keys);

    MPI_Finalize();
}




//spatial decomp with load balancing
//build local tree, reduce, broadcast
void NBodyBH4_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Alloc data
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    long currAlloc = 2*localN;

    int startN[world_size];
    int start3N[world_size];
    int numN[world_size];
    int num3N[world_size];

    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    //for local proc's bodies
    err = allocData_NB(currAlloc, &r, &v, &a, &m, &work);
    if (err) {
        fprintf(stderr, "Could not alloc local data for N=%ld\n", localN);
        exit(ALLOC_ERROR);
    }


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * Init data
 **********************************/

    spatialKey_t* keys = NULL;
    long* idx = NULL, *idx2 = NULL;
    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL, *wwork = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData_NB(N, &rr, &vv, &aa, &mm, &wwork);
        err |= allocSpatialKeys_NB(N, &keys);
        //if malloc returns null, set err bit
        err |= !(idx = (long*) malloc(sizeof(long)*N));
        err |= !(idx2 = (long*) malloc(sizeof(long)*N));

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);
        for (long i = 0; i < N; ++i) {
            wwork[i] = N;
        }

        //sort before scattering
        computeSpatialKeys_NB(N, rr, keys);
        sortSpatialKeys_NB(N, keys, idx);
#if NBODY_SIM_WITH_RENDERER
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap4N_NB(N, idx2, colors);
#endif
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, idx2, rr);
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap3N_NB(N, idx2, vv);
        memcpy(idx2, idx, sizeof(long)*N);
        sortByIdxMap_NB(N, idx2, mm);
        sortByIdxMap_NB(N, idx, wwork);

        computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }

    err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size; ++i) {
        num3N[i] = 3*numN[i];
        start3N[i] = 3*startN[i];
    }

    if (numN[world_rank] > currAlloc) {
        while (numN[world_rank] > currAlloc) {
            currAlloc += (N / world_size);
        }
        err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
    }

    localN = numN[world_rank];
    err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                        r, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                        v, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(aa, num3N, start3N, MPI_DOUBLE,
                        a, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                        m, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(wwork, numN, startN, MPI_DOUBLE,
                        work, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        //gather for sorting and work partition. We don't need accels so don't gather that.
        err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                          rr, num3N, start3N, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
        err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                          vv, num3N, start3N, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
        err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                          mm, numN, startN, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
        err = MPI_Gatherv( work, localN, MPI_DOUBLE,
                          wwork, numN, startN, MPI_DOUBLE,
                          0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            //compute keys, sort them, then sort each data array by index map
            //this avoids simultaneous sorting SoA's or having a AoS.
            computeSpatialKeys_NB(N, rr, keys);
            sortSpatialKeys_NB(N, keys, idx);
#if NBODY_SIM_WITH_RENDERER
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap4N_NB(N, idx2, colors);
#endif
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, rr);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, vv);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap_NB(N, idx2, mm);
            sortByIdxMap_NB(N, idx, wwork);

            computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);
        }

        err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
        err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 0; i < world_size; ++i) {
            num3N[i] = 3*numN[i];
            start3N[i] = 3*startN[i];
        }

        if (numN[world_rank] > currAlloc) {
            while (numN[world_rank] > currAlloc) {
                currAlloc += (N / world_size);
            }
            err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
        }

        localN = numN[world_rank];
        err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                            r, num3N[world_rank], MPI_DOUBLE,
                            0, MPI_COMM_WORLD);
        err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                            v, num3N[world_rank], MPI_DOUBLE,
                            0, MPI_COMM_WORLD);
        err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                            m, numN[world_rank], MPI_DOUBLE,
                            0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            domainSize = computeDomainSize_NB(N, rr);
        }
        MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        mapReduceBuildOctrees_NBMPI(r, m, localN, domainSize, &octree, &treeData, &treeAlloc);

        if (world_rank == 0) {
            treeSize = serializeOctree_NB(octree, &treeData, &treeAlloc);
        }

        MPI_Bcast(&treeSize, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (treeSize > treeAlloc) {
                treeData = (uint8_t*) realloc(treeData, sizeof(uint8_t)*treeSize);
                treeAlloc = treeSize;
            }
        }
        MPI_Bcast(treeData, treeSize, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            freeOctree_NB(octree);
            octree = deserializeOctree_NB(treeData, treeSize);
        }

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesOctreeBH_NB(localN, m, r, work, a, octree, tmp1_node, tmp2_node, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                      rr, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                      vv, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                      mm, numN, startN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(work);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(wwork);
    free(tmp1_node);
    free(tmp2_node);
    free(idx);
    free(idx2);
    free(keys);

    MPI_Finalize();
}




//distributed spatial decomp and load balancing
//build local tree, reduce, broadcast
void NBodyBH5_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init MPI
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    long currAlloc = 2*localN;

    int startN[world_size];
    int start3N[world_size];
    int numN[world_size];
    int num3N[world_size];
    int minN = localN;


/**********************************
 * Allocate processor's local data
 **********************************/

    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    long *localIdx1 = NULL, *localIdx2 = NULL;
    spatialKey_t* keys = NULL;
    err = allocData_NB(currAlloc, &r, &v, &a, &m, &work);
    err |= allocSpatialKeys_NB(currAlloc, &keys);
    err |= !(localIdx1 = (long*) malloc(sizeof(long)*currAlloc));
    err |= !(localIdx2 = (long*) malloc(sizeof(long)*currAlloc));

    if (err) {
        fprintf(stderr, "Could not alloc local data for N=%ld\n", localN);
        exit(ALLOC_ERROR);
    }


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * On the root node, initialize all the starting data.
 **********************************/

    spatialKey_t* kkeys = NULL;
    long* idx = NULL, *idx2 = NULL;
    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL, *wwork = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData_NB(N, &rr, &vv, &aa, &mm, &wwork);

        // if (1) {
        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
            err |= allocSpatialKeys_NB(N, &kkeys);
            //if malloc returns null, set err bit
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#if NBODY_SIM_WITH_RENDERER
        else {
            err |= allocSpatialKeys_NB(N, &kkeys);
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#endif

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);
        for (long i = 0; i < N; ++i) {
            wwork[i] = N;
        }


#if NBODY_SIM_WITH_RENDERER
        {
#else
        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
#endif
            //sort before scattering
            computeSpatialKeys_NB(N, rr, kkeys);
            sortSpatialKeys_NB(N, kkeys, idx);
            memcpy(idx2, idx, sizeof(long)*N);
#if NBODY_SIM_WITH_RENDERER
            sortByIdxMap4N_NB(N, idx2, colors);
            memcpy(idx2, idx, sizeof(long)*N);
#endif
            sortByIdxMap3N_NB(N, idx2, rr);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, vv);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap_NB(N, idx2, mm);
            sortByIdxMap_NB(N, idx, wwork);
        }

        computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }


/**********************************
 * Scatter the root-initialized data to all processors.
 **********************************/
    err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size; ++i) {
        num3N[i] = 3*numN[i];
        start3N[i] = 3*startN[i];
    }

    if (numN[world_rank] > currAlloc) {
        while (numN[world_rank] > currAlloc) {
            currAlloc += (N / world_size);
        }
        err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
        err |= reallocSpatialKeys_NB(currAlloc, &keys);
        err |= !(localIdx1 = (long*) realloc(localIdx1, sizeof(long)*currAlloc));
        err |= !(localIdx2 = (long*) realloc(localIdx2, sizeof(long)*currAlloc));
    }

    localN = numN[world_rank];
    err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                        r, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                        v, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(aa, num3N, start3N, MPI_DOUBLE,
                        a, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                        m, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(wwork, numN, startN, MPI_DOUBLE,
                        work, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

    double domainSize;
    domainSize = computeDomainSize_NB(localN, r);
    MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (N >=  2*SORT_EXCHANGE_SIZE*world_size) {
        //If large enough, haven't yet sorted, do so after scattering
        computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);
    }


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        // if (1) {
        if (minN < 2*SORT_EXCHANGE_SIZE) {
            if (kkeys == NULL) {
                err |= allocSpatialKeys_NB(N, &kkeys);
            }
            if (idx == NULL) {
                //if malloc returns null, set err bit
                err |= !(idx = (long*) malloc(sizeof(long)*N));
            }
            if (idx2 == NULL) {
                err |= !(idx2 = (long*) malloc(sizeof(long)*N));
            }
            if (err) {
                fprintf(stderr, "Could not alloc data for N=%ld\n", N);
                exit(ALLOC_ERROR);
            }

            //gather for sorting and work partition. We don't need accels so don't gather that.
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                              rr, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                              vv, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                              mm, numN, startN, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( work, localN, MPI_DOUBLE,
                              wwork, numN, startN, MPI_DOUBLE,
                              0, MPI_COMM_WORLD);

            if (world_rank == 0) {

                //compute keys, sort them, then sort each data array by index map
                //this avoids simultaneous sorting SoA's or having a AoS.
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
    #if NBODY_SIM_WITH_RENDERER
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap4N_NB(N, idx2, colors);
    #endif
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, vv);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap_NB(N, idx2, mm);
                sortByIdxMap_NB(N, idx, wwork);

                computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 0; i < world_size; ++i) {
                num3N[i] = 3*numN[i];
                start3N[i] = 3*startN[i];
            }

            if (numN[world_rank] > currAlloc) {
                while (numN[world_rank] > currAlloc) {
                    currAlloc += (N / world_size);
                }
                err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
                err |= reallocSpatialKeys_NB(currAlloc, &keys);
            }

            localN = numN[world_rank];
            err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                                r, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                                v, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                                m, numN[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                domainSize = computeDomainSize_NB(N, rr);
            }
            MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
#if NBODY_SIM_WITH_RENDERER
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                          rr, num3N, start3N, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                sortByIdxMap4N_NB(N, idx, colors);
            }
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            //PARALLEL COMPUTATION OF OCTREE EXTENTS.
            domainSize = computeDomainSize_NB(localN, r);
            MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            //PARALLEL SORT
            computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);

            //LOAD BALANCE
            distributeWork_NBMPI(N, &localN, &currAlloc, &r, &v, &a, &m, &work, &keys, &localIdx1, &localIdx2);

            //Share size of every proc with every proc for future Gatherv.
            MPI_Allgather(&localN, 1, MPI_INT,
                          numN, 1, MPI_INT, MPI_COMM_WORLD);
            num3N[0] = numN[0]*3;
            startN[0] = 0;
            start3N[0] = 0;
            minN = numN[0];
            for (int i = 1; i < world_size; ++i) {
                num3N[i] = numN[i]*3;
                startN[i] = startN[i-1] + numN[i-1];
                start3N[i] = start3N[i-1] + num3N[i-1];
                minN = numN[i] < minN ? numN[i] : minN;
            }
        }


        mapReduceBuildOctrees_NBMPI(r, m, localN, domainSize, &octree, &treeData, &treeAlloc);

        if (world_rank == 0) {
            treeSize = serializeOctree_NB(octree, &treeData, &treeAlloc);
        }

        MPI_Bcast(&treeSize, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (treeSize > treeAlloc) {
                treeData = (uint8_t*) realloc(treeData, sizeof(uint8_t)*treeSize);
                treeAlloc = treeSize;
            }
        }

        MPI_Bcast(treeData, treeSize, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            freeOctree_NB(octree);
            octree = deserializeOctree_NB(treeData, treeSize);
        }

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesOctreeBH_NB(localN, m, r, work, a, octree, tmp1_node, tmp2_node, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);
    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Allgather(&localN, 1, MPI_INT,
                  numN, 1, MPI_INT, MPI_COMM_WORLD);
    num3N[0] = numN[0]*3;
    startN[0] = 0;
    start3N[0] = 0;
    for (int i = 1; i < world_size; ++i) {
        num3N[i] = numN[i]*3;
        startN[i] = startN[i-1] + numN[i-1];
        start3N[i] = start3N[i-1] + num3N[i-1];
    }

    err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                      rr, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                      vv, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                      mm, numN, startN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(work);
    free(keys);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(wwork);
    free(kkeys);
    free(tmp1_node);
    free(tmp2_node);
    free(localIdx1);
    free(localIdx2);
    free(idx);
    free(idx2);

    MPI_Finalize();
}


//distributed spatial decomp and load balancing
//build local tree, broadcast branch nodes,
//dynamic request data during traversal
void NBodyBH6_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init MPI
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    long currAlloc = 2*localN;

    int startN[world_size];
    int start3N[world_size];
    int numN[world_size];
    int num3N[world_size];
    int minN = localN;


/**********************************
 * Allocate processor's local data
 **********************************/

    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    long *localIdx1 = NULL, *localIdx2 = NULL;
    spatialKey_t* keys = NULL;
    err = allocData_NB(currAlloc, &r, &v, &a, &m, &work);
    err |= allocSpatialKeys_NB(currAlloc, &keys);
    err |= !(localIdx1 = (long*) malloc(sizeof(long)*currAlloc));
    err |= !(localIdx2 = (long*) malloc(sizeof(long)*currAlloc));

    if (err) {
        fprintf(stderr, "Could not alloc local data for N=%ld\n", localN);
        exit(ALLOC_ERROR);
    }


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * On the root node, initialize all the starting data.
 **********************************/

    spatialKey_t* kkeys = NULL;
    long* idx = NULL, *idx2 = NULL;
    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL, *wwork = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData_NB(N, &rr, &vv, &aa, &mm, &wwork);

        // if (1) {
        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
            err |= allocSpatialKeys_NB(N, &kkeys);
            //if malloc returns null, set err bit
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#if NBODY_SIM_WITH_RENDERER
        else {
            err |= allocSpatialKeys_NB(N, &kkeys);
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#endif

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);
        for (long i = 0; i < N; ++i) {
            wwork[i] = N;
        }

        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
            //sort before scattering
            computeSpatialKeys_NB(N, rr, kkeys);
            sortSpatialKeys_NB(N, kkeys, idx);
#if NBODY_SIM_WITH_RENDERER
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap4N_NB(N, idx2, colors);
#endif
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, rr);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, vv);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap_NB(N, idx2, mm);
            sortByIdxMap_NB(N, idx, wwork);
        }

        computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }


/**********************************
 * Scatter the root-initialized data to all processors.
 **********************************/
    err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size; ++i) {
        num3N[i] = 3*numN[i];
        start3N[i] = 3*startN[i];
    }

    if (numN[world_rank] > currAlloc) {
        while (numN[world_rank] > currAlloc) {
            currAlloc += (N / world_size);
        }
        err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
        err |= reallocSpatialKeys_NB(currAlloc, &keys);
        err |= !(localIdx1 = (long*) realloc(localIdx1, sizeof(long)*currAlloc));
        err |= !(localIdx2 = (long*) realloc(localIdx2, sizeof(long)*currAlloc));
    }

    localN = numN[world_rank];
    err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                        r, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                        v, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(aa, num3N, start3N, MPI_DOUBLE,
                        a, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                        m, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(wwork, numN, startN, MPI_DOUBLE,
                        work, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

    double domainSize;
    domainSize = computeDomainSize_NB(localN, r);
    MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (N >= 2*SORT_EXCHANGE_SIZE*world_size) {
        computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);
    }


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    NBodyHOT_t HOT = {NULL, 0, 0};
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        // if (1) {
        if (minN < 2*SORT_EXCHANGE_SIZE) {
            if (kkeys == NULL) {
                err |= allocSpatialKeys_NB(N, &kkeys);
            }
            if (idx == NULL) {
                //if malloc returns null, set err bit
                err |= !(idx = (long*) malloc(sizeof(long)*N));
            }
            if (idx2 == NULL) {
                err |= !(idx2 = (long*) malloc(sizeof(long)*N));
            }
            if (err) {
                fprintf(stderr, "Could not alloc data for N=%ld\n", N);
                exit(ALLOC_ERROR);
            }

            //gather for sorting and work partition. We don't need accels so don't gather that.
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                              rr, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                              vv, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                              mm, numN, startN, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( work, localN, MPI_DOUBLE,
                              wwork, numN, startN, MPI_DOUBLE,
                              0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                //compute keys, sort them, then sort each data array by index map
                //this avoids simultaneous sorting SoA's or having a AoS.
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
    #if NBODY_SIM_WITH_RENDERER
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap4N_NB(N, idx2, colors);
    #endif
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, vv);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap_NB(N, idx2, mm);
                sortByIdxMap_NB(N, idx, wwork);

                computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);
            }

            err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 0; i < world_size; ++i) {
                num3N[i] = 3*numN[i];
                start3N[i] = 3*startN[i];
            }

            if (numN[world_rank] > currAlloc) {
                while (numN[world_rank] > currAlloc) {
                    currAlloc += (N / world_size);
                }
                err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
                err |= reallocSpatialKeys_NB(currAlloc, &keys);
            }

            localN = numN[world_rank];
            err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                                r, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                                v, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                                m, numN[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                domainSize = computeDomainSize_NB(N, rr);
            }
            MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
#if NBODY_SIM_WITH_RENDERER
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                          rr, num3N, start3N, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                sortByIdxMap4N_NB(N, idx, colors);
            }
#endif
            //PARALLEL COMPUTATION OF OCTREE EXTENTS.
            domainSize = computeDomainSize_NB(localN, r);
            MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            //PARALLEL SORT
            computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);

            //LOAD BALANCE
            distributeWork_NBMPI(N, &localN, &currAlloc, &r, &v, &a, &m, &work, &keys, &localIdx1, &localIdx2);

            //Share size of every proc with every proc for future Gatherv.
            MPI_Allgather(&localN, 1, MPI_INT,
                          numN, 1, MPI_INT, MPI_COMM_WORLD);
            num3N[0] = numN[0]*3;
            startN[0] = 0;
            start3N[0] = 0;
            minN = numN[0];
            for (int i = 1; i < world_size; ++i) {
                num3N[i] = numN[i]*3;
                startN[i] = startN[i-1] + numN[i-1];
                start3N[i] = start3N[i-1] + num3N[i-1];
                minN = numN[i] < minN ? numN[i] : minN;
            }
        }


        HOT = buildHashedOctreeInPlace_NB(localN, r, m, domainSize, HOT);

        exchangeBorderBodies_NB(HOT, localN, r, m);

        exchangeBranchNodes_NB(HOT, &treeData, &treeAlloc);

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            //Note that masses doesn't change here because we aren't doing any sorting.
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesDistributedHOTBH_NB(localN, m, r, work, a, HOT,
            (NBodyHOTNode_t**) tmp1_node, (NBodyHOTNode_t**) tmp2_node, &treeData, &treeAlloc, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);

    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Allgather(&localN, 1, MPI_INT,
                  numN, 1, MPI_INT, MPI_COMM_WORLD);
    num3N[0] = numN[0]*3;
    startN[0] = 0;
    start3N[0] = 0;
    for (int i = 1; i < world_size; ++i) {
        num3N[i] = numN[i]*3;
        startN[i] = startN[i-1] + numN[i-1];
        start3N[i] = start3N[i-1] + num3N[i-1];
    }

    err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                      rr, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                      vv, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                      mm, numN, startN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(work);
    free(keys);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(wwork);
    free(kkeys);
    free(tmp1_node);
    free(tmp2_node);
    free(localIdx1);
    free(localIdx2);
    free(idx);
    free(idx2);

    MPI_Finalize();
}



//distributed spatial decomp and load balancing
//build local tree, broadcast branch nodes,
//dynamic and asynchronous request data during traversal
void NBodyBH7_NBMPI(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init MPI
 **********************************/

    MPI_Init(NULL, NULL);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int err;
    long localN = N / world_size; //assume divisible
    long currAlloc = 2*localN;

    int startN[world_size];
    int start3N[world_size];
    int numN[world_size];
    int num3N[world_size];
    int minN = localN;


/**********************************
 * Allocate processor's local data
 **********************************/

    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    long *localIdx1 = NULL, *localIdx2 = NULL;
    spatialKey_t* keys = NULL;
    err = allocData_NB(currAlloc, &r, &v, &a, &m, &work);
    err |= allocSpatialKeys_NB(currAlloc, &keys);
    err |= !(localIdx1 = (long*) malloc(sizeof(long)*currAlloc));
    err |= !(localIdx2 = (long*) malloc(sizeof(long)*currAlloc));

    if (err) {
        fprintf(stderr, "Could not alloc local data for N=%ld\n", localN);
        exit(ALLOC_ERROR);
    }


/**********************************
 * Start renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = NULL;
    NBodyRenderer* renderer = NULL;
    if (world_rank == 0) {
        colors = createColors(N);
        renderer = new NBodyRenderer();
        renderer->startRenderThread();
        // renderer.updatePoints(N, r);

        // renderer.setOctree(octree);
        // renderer.joinRenderThread();
        // return 0;
    }
#endif


/**********************************
 * On the root node, initialize all the starting data.
 **********************************/

    spatialKey_t* kkeys = NULL;
    long* idx = NULL, *idx2 = NULL;
    double *rr = NULL, *vv = NULL, *aa = NULL, *mm = NULL, *wwork = NULL;
    double Epot, Ekin, E0;
    if (world_rank == 0) {
        //init total data as root

        err = allocData_NB(N, &rr, &vv, &aa, &mm, &wwork);

        // if (1) {
        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
            err |= allocSpatialKeys_NB(N, &kkeys);
            //if malloc returns null, set err bit
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#if NBODY_SIM_WITH_RENDERER
        else {
            err |= allocSpatialKeys_NB(N, &kkeys);
            err |= !(idx = (long*) malloc(sizeof(long)*N));
            err |= !(idx2 = (long*) malloc(sizeof(long)*N));
        }
#endif

        if (err) {
            fprintf(stderr, "Could not alloc data for N=%ld\n", N);
            exit(ALLOC_ERROR);
        }

        err = initData_NB(N, seed, rr, vv, aa, mm);
        for (long i = 0; i < N; ++i) {
            wwork[i] = N;
        }

        if (N < 2*SORT_EXCHANGE_SIZE*world_size) {
            //sort before scattering
            computeSpatialKeys_NB(N, rr, kkeys);
            sortSpatialKeys_NB(N, kkeys, idx);
#if NBODY_SIM_WITH_RENDERER
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap4N_NB(N, idx2, colors);
#endif
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, rr);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap3N_NB(N, idx2, vv);
            memcpy(idx2, idx, sizeof(long)*N);
            sortByIdxMap_NB(N, idx2, mm);
            sortByIdxMap_NB(N, idx, wwork);
        }

        computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);

        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "E0: %.15g\n", E0);
    }


/**********************************
 * Scatter the root-initialized data to all processors.
 **********************************/
    err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size; ++i) {
        num3N[i] = 3*numN[i];
        start3N[i] = 3*startN[i];
    }

    if (numN[world_rank] > currAlloc) {
        while (numN[world_rank] > currAlloc) {
            currAlloc += (N / world_size);
        }
        err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
        err |= reallocSpatialKeys_NB(currAlloc, &keys);
        err |= !(localIdx1 = (long*) realloc(localIdx1, sizeof(long)*currAlloc));
        err |= !(localIdx2 = (long*) realloc(localIdx2, sizeof(long)*currAlloc));
    }

    localN = numN[world_rank];
    err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                        r, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                        v, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(aa, num3N, start3N, MPI_DOUBLE,
                        a, num3N[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                        m, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
    err = MPI_Scatterv(wwork, numN, startN, MPI_DOUBLE,
                        work, numN[world_rank], MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

    double domainSize;
    domainSize = computeDomainSize_NB(localN, r);
    MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (N >= 2*SORT_EXCHANGE_SIZE*world_size) {
        computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);
    }


/**********************************
 * Start sim
 **********************************/
    double dt_out = 0.1;
    double t_out = dt_out;

    NBodyHOT_t HOT = {NULL, 0, 0};
    NBOctree_t* octree = NULL;
    int updateInplace = 0;
    uint8_t* treeData = NULL;
    size_t treeSize = 0, treeAlloc = 0;

    //working space needs to be N, not local N, because that is max number of interactions
    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);

    unsigned long long startTimer;
    _startTimerParallel(&startTimer);
    for (double t = 0.0; t < t_end; t += dt) {
// #if NBODY_SIM_WITH_RENDERER
        // if (renderer.shouldClose()) {
            // break;
        // }
// #endif

        //kick and drift
        performNBodyHalfStepA(localN, dt, r, v, a, m);

        // if (1) {
        if (minN < 2*SORT_EXCHANGE_SIZE) {
            if (kkeys == NULL) {
                err |= allocSpatialKeys_NB(N, &kkeys);
            }
            if (idx == NULL) {
                //if malloc returns null, set err bit
                err |= !(idx = (long*) malloc(sizeof(long)*N));
            }
            if (idx2 == NULL) {
                err |= !(idx2 = (long*) malloc(sizeof(long)*N));
            }
            if (err) {
                fprintf(stderr, "Could not alloc data for N=%ld\n", N);
                exit(ALLOC_ERROR);
            }

            //gather for sorting and work partition. We don't need accels so don't gather that.
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                              rr, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                              vv, num3N, start3N, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                              mm, numN, startN, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
            err = MPI_Gatherv( work, localN, MPI_DOUBLE,
                              wwork, numN, startN, MPI_DOUBLE,
                              0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                //compute keys, sort them, then sort each data array by index map
                //this avoids simultaneous sorting SoA's or having a AoS.
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
    #if NBODY_SIM_WITH_RENDERER
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap4N_NB(N, idx2, colors);
    #endif
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, vv);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap_NB(N, idx2, mm);
                sortByIdxMap_NB(N, idx, wwork);

                computeWorkPartitions_NBMPI(N, wwork, startN, numN, world_size);
            }

            err = MPI_Bcast(numN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            err = MPI_Bcast(startN, world_size, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 0; i < world_size; ++i) {
                num3N[i] = 3*numN[i];
                start3N[i] = 3*startN[i];
            }

            if (numN[world_rank] > currAlloc) {
                while (numN[world_rank] > currAlloc) {
                    currAlloc += (N / world_size);
                }
                err = reallocData_NB(currAlloc, &r, &v, &a, &m, &work);
                err |= reallocSpatialKeys_NB(currAlloc, &keys);
            }

            localN = numN[world_rank];
            err = MPI_Scatterv(rr, num3N, start3N, MPI_DOUBLE,
                                r, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(vv, num3N, start3N, MPI_DOUBLE,
                                v, num3N[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);
            err = MPI_Scatterv(mm, numN, startN, MPI_DOUBLE,
                                m, numN[world_rank], MPI_DOUBLE,
                                0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                domainSize = computeDomainSize_NB(N, rr);
            }
            MPI_Bcast(&domainSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {

#if NBODY_SIM_WITH_RENDERER
            err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                          rr, num3N, start3N, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);

            if (world_rank == 0) {
                computeSpatialKeys_NB(N, rr, kkeys);
                sortSpatialKeys_NB(N, kkeys, idx);
                memcpy(idx2, idx, sizeof(long)*N);
                sortByIdxMap3N_NB(N, idx2, rr);
                sortByIdxMap4N_NB(N, idx, colors);
            }
#endif
            //PARALLEL COMPUTATION OF OCTREE EXTENTS.
            domainSize = computeDomainSize_NB(localN, r);
            MPI_Allreduce(MPI_IN_PLACE, &domainSize, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            //PARALLEL SORT
            computeAndSortKeys_NBMPI(N, localN, domainSize, keys, r, v, m, work, localIdx1, localIdx2);

            //LOAD BALANCE
            distributeWork_NBMPI(N, &localN, &currAlloc, &r, &v, &a, &m, &work, &keys, &localIdx1, &localIdx2);

            //Share size of every proc with every proc for future Gatherv.
            MPI_Allgather(&localN, 1, MPI_INT,
                          numN, 1, MPI_INT, MPI_COMM_WORLD);
            num3N[0] = numN[0]*3;
            startN[0] = 0;
            start3N[0] = 0;
            minN = numN[0];
            for (int i = 1; i < world_size; ++i) {
                num3N[i] = numN[i]*3;
                startN[i] = startN[i-1] + numN[i-1];
                start3N[i] = start3N[i-1] + num3N[i-1];
                minN = numN[i] < minN ? numN[i] : minN;
            }
        }


        HOT = buildHashedOctreeInPlace_NB(localN, r, m, domainSize, HOT);

        exchangeBorderBodies_NB(HOT, localN, r, m);

        exchangeBranchNodes_NB(HOT, &treeData, &treeAlloc);

#if NBODY_SIM_WITH_RENDERER
        if (world_rank == 0) {
            //update renderer before next kick since positions only change in step A
            //Note that masses doesn't change here because we aren't doing any sorting.
            if (renderer->needsUpdate()) {
                renderer->updatePoints(N, rr, mm, 1, colors);
            }
        }
#endif

        //update accels
        computeForcesAsyncDistributedHOTBH_NB(localN, m, r, work, a, HOT,
            (NBodyHOTNode_t**) tmp1_node, (NBodyHOTNode_t**) tmp2_node, &treeData, &treeAlloc, theta);

        //second kick
        performNBodyHalfStepB(localN, dt, r, v, a, m);

    }

    float elapsed = 0.0;
    _stopTimerAddElapsedParallel(&startTimer, &elapsed);

    MPI_Allgather(&localN, 1, MPI_INT,
                  numN, 1, MPI_INT, MPI_COMM_WORLD);
    num3N[0] = numN[0]*3;
    startN[0] = 0;
    start3N[0] = 0;
    for (int i = 1; i < world_size; ++i) {
        num3N[i] = numN[i]*3;
        startN[i] = startN[i-1] + numN[i-1];
        start3N[i] = start3N[i-1] + num3N[i-1];
    }

    err = MPI_Gatherv( r, localN*3, MPI_DOUBLE,
                      rr, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( v, localN*3, MPI_DOUBLE,
                      vv, num3N, start3N, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);
    err = MPI_Gatherv( m, localN, MPI_DOUBLE,
                      mm, numN, startN, MPI_DOUBLE,
                       0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        Epot = computeEpot_NB(N, mm, rr);
        Ekin = computeEkin_NB(N, mm, vv);
        E0 = Epot + Ekin;
        fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
        fprintf(stderr, "Eend: %.15g\n", E0);
        fprintf(stderr, "Elapsed Time: %.15g\n", elapsed);
    }

#if NBODY_SIM_WITH_RENDERER
    if (world_rank == 0) {
        //End of simulation, wait for renderer to finish.
        renderer->joinRenderThread();
        //End of simultation, force renderer to finish.
        // renderer.stopRenderThread();
        delete renderer;
    }
#endif

/**********************************
 * Clean up
 **********************************/
    free(r);
    free(v);
    free(a);
    free(m);
    free(work);
    free(keys);
    free(rr);
    free(vv);
    free(aa);
    free(mm);
    free(wwork);
    free(kkeys);
    free(tmp1_node);
    free(tmp2_node);
    free(localIdx1);
    free(localIdx2);
    free(idx);
    free(idx2);

    MPI_Finalize();
}
