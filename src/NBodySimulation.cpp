
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

#if NBODY_PARALLEL
#include "NBodyParallel.h"
int ExecutorThreadPool::maxThreads = NBODY_NPROCS;
#endif

#include "Unix_Timer.h"



void NBodySimSerial(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init data
 **********************************/
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    int err = allocData(N, &r, &v, &a, &m, &work);
    if (err) {
        fprintf(stderr, "Could not alloc data for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }

    // N = 3;
    // err = initData3BodyFigureEight(r, v, a, m);
    err = initData(N, seed, r, v, a, m);

    //for first round just guess work = N;
    for (long i = 0; i < N; ++i) {
        work[i] = N;
    }

    double Epot = computeEpot_NB(N, m, r);
    double Ekin = computeEkin_NB(N, m, v);
    double E0 = Epot + Ekin;
    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "E0: %.15g\n", E0);


/**********************************
 * Temporary working space
 **********************************/
    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace;

    const NBOctreeNode_t** tmp1_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    const NBOctreeNode_t** tmp2_node = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    if (tmp1_node == NULL || tmp2_node == NULL) {
        fprintf(stderr, "Could not alloc enough working space for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }


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

        //build Octree and compute forces to update acceleration.
        domainSize = computeDomainSize_NB(N, r);
        updateInplace = buildOctreeInPlace_NB(N, r, m, domainSize, &octree);

#if NBODY_SIM_WITH_RENDERER
        //update renderer before next kick since positions only change in step A
        if (renderer.needsUpdate()) {
            renderer.updatePoints(N, r, m, domainSize, colors);
        }
#endif

        //update accels
        // computeForces(N, m, r, a);
        // computeForcesMonoOctreeBH_NB(N, m, r, work, a, octree, tmp1_node, tmp2_node, theta);
        computeForcesOctreeBH_NB(N, m, r, work, a, octree, tmp1_node, tmp2_node, theta);

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
    freeData(r, v, a, m, work);
    free(tmp1_node);
    free(tmp2_node);
}



#if NBODY_PARALLEL
void NBodySimParallel(long N, double dt, double t_end, time_t seed, double theta) {

/**********************************
 * Init data
 **********************************/
    double *r = NULL, *v = NULL, *a = NULL, *m = NULL, *work = NULL;
    int err = allocData(N, &r, &v, &a, &m, &work);
    if (err) {
        fprintf(stderr, "Could not alloc data for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }

    // N = 3;
    // err = initData3BodyFigureEight(r, v, a, m);
    err = initData(N, seed, r, v, a, m);

    //for first round just guess work = N;
    for (long i = 0; i < N; ++i) {
        work[i] = N;
    }
    long startN[NBODY_NPROCS];
    long numN[NBODY_NPROCS];
    computeWorkPartitions(N, work, startN, numN, NBODY_NPROCS);

    double Epot = computeEpot_NB(N, m, r);
    double Ekin = computeEkin_NB(N, m, v);
    double E0 = Epot + Ekin;
    fprintf(stderr, "Ekin: %.15g\nEpot: %.15g\n", Ekin, Epot);
    fprintf(stderr, "E0: %.15g\n", E0);


/**********************************
 * Allocate temporary working space
 **********************************/
    spatialKey_t* keys;
    err = allocSpatialKeys_NB(N, &keys);
    if (err) {
        fprintf(stderr, "Could not alloc keys for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }

    double domainSize;
    NBOctree_t* octree = NULL;
    int updateInplace;

    long* idx = (long*) malloc(sizeof(long)*N);
    //pointers and longs are same size so let's hack our temporary working space
    const NBOctreeNode_t** tmpNodes1[NBODY_NPROCS];
    const NBOctreeNode_t** tmpNodes2[NBODY_NPROCS];
    for (int i = 0; i < NBODY_NPROCS; ++i) {
        //interleaving like this actually causes better locality
        //where thread i uses both tmpNodes1[i] and tmpNodes2[i]
        tmpNodes1[i] = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
        tmpNodes2[i] = (const NBOctreeNode_t**) malloc(sizeof(NBOctreeNode_t*)*N);
    }

    NBOctree_t* trees[NBODY_NPROCS];
    for (int i = 0; i < NBODY_NPROCS; ++i) {
        trees[i] = NULL;
    }

    long* tmpLong[NBODY_NPROCS];
    if (sizeof(long) == sizeof(NBOctreeNode_t**)) {
        for (int i = 0; i < NBODY_NPROCS; ++i) {
            tmpLong[i] = (long*) tmpNodes1[i];
        }
    } else {
        for (int i = 0; i < NBODY_NPROCS; ++i) {
            tmpLong[i] = (long*) malloc(sizeof(long)*N);
        }
    }

    if (idx == NULL ||
        tmpNodes1[NBODY_NPROCS-1] == NULL ||
        tmpNodes2[NBODY_NPROCS-1] == NULL ||
        tmpLong[NBODY_NPROCS-1] == NULL)
    {
        fprintf(stderr, "Could not alloc enough working space for N=%ld\n", N);
        exit(ALLOC_ERROR);
    }


/**********************************
 * Setup renderer
 **********************************/
#if NBODY_SIM_WITH_RENDERER
    float* colors = createColors(N);
    NBodyRenderer renderer;
    renderer.startRenderThread();
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
        //kick, drift
        performNBodyHalfStepAParallel_NB(dt, r, v, a, m, startN, numN);

        //get spatial ordering, sort, and then partition.
        //not worth parallelizing key gen, only 0.002s for 100k keys.
        computeSpatialKeys_NB(N, r, keys);
        sortSpatialKeys_NB(N, keys, idx);
        parallelSortByIdxMap_NB(
            N, r, v, work,
        #if NBODY_SIM_WITH_RENDERER
            colors,
        #endif
            idx,
            tmpLong
        );
        computeWorkPartitions(N, work, startN, numN, NBODY_NPROCS);

        //build Octree
        domainSize = computeDomainSize_NB(N, r);
        updateInplace = mapReduceBuildOctreesInPlace_NB(r, m, domainSize, trees, startN, numN);
        octree = trees[0];

#if NBODY_SIM_WITH_RENDERER
        //update renderer before next kick since positions only change in step A
        if (renderer.needsUpdate()) {
            renderer.updatePoints(N, r, m, domainSize, colors);
        }
#endif
        //compute forces to update acceleration.
        computeForcesOctreeBHParallel_NB(m, r, work, a, octree, tmpNodes1, tmpNodes2, theta, startN, numN);

        //kick
        performNBodyHalfStepBParallel_NB(dt, r, v, a, m, startN, numN);
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
    freeData(r, v, a, m, work);
    freeSpatialKeys_NB(keys);
    free(idx);
    for (int i = 0; i < NBODY_NPROCS; ++i) {
        free(tmpNodes1[i]);
        free(tmpNodes2[i]);
    }
    if (sizeof(long) != sizeof(NBOctreeNode_t**)) {
        for (int i = 0; i < NBODY_NPROCS; ++i) {
            free(tmpLong[i]);
        }
    }

}
#endif
