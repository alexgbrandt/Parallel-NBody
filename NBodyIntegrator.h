
#ifndef _NBODY_INTEGRATOR_
#define _NBODY_INTEGRATOR_

#ifdef __cplusplus
extern "C" {
#endif


#include "NBodyForces.h"

/**
 * Perform the first half of the leafprog integration,
 * "kick, drift". Updating velocities to the half step
 * and positions to the full step.
 *
 * @param n, the number of bodies.
 * @param dt, the time step for integration.
 * @param r, an array of 3*n doubles holding the positions of the bodies.
 * @param v, an array of 3*n doubles holding the velocitites of the bodies.
 * @param a, an array of 3*n doubles holding the acceleration of the bodies.
 * @param m, an array of n doubles holding the masses of the bodies.
 */
static inline void performNBodyHalfStepA(long n, double dt,
    double* __restrict__ r,
    double* __restrict__ v,
    const double* __restrict__ a,
    const double* __restrict__ m)
{
    for (long i = 0; i < n; ++i) {
        //kick, drift
        v[3*i + 0] += 0.5 * a[3*i + 0] * dt;
        r[3*i + 0] += v[3*i + 0] * dt;
        v[3*i + 1] += 0.5 * a[3*i + 1] * dt;
        r[3*i + 1] += v[3*i + 1] * dt;
        v[3*i + 2] += 0.5 * a[3*i + 2] * dt;
        r[3*i + 2] += v[3*i + 2] * dt;
    }
}


/**
 * Perform the second half of the leafprog integration,
 * "kick2". Updating velocities from the half step to
 * the full step.
 * Positions are not actually updated here but included
 * for symmetry.
 *
 * @param n, the number of bodies.
 * @param dt, the time step for integration.
 * @param r, an array of 3*n doubles holding the positions of the bodies.
 * @param v, an array of 3*n doubles holding the velocitites of the bodies.
 * @param a, an array of 3*n doubles holding the acceleration of the bodies.
 * @param m, an array of n doubles holding the masses of the bodies.
 */
static inline void performNBodyHalfStepB(long n, double dt,
    const double* __restrict__ r,
    double* __restrict__ v,
    const double* __restrict__ a,
    const double* __restrict__ m)
{
    for (long i = 0; i < n; ++i) {
        //kick
        v[3*i + 0] += 0.5 * a[3*i + 0] * dt;
        v[3*i + 1] += 0.5 * a[3*i + 1] * dt;
        v[3*i + 2] += 0.5 * a[3*i + 2] * dt;
    }
}



#ifdef __cplusplus
}
#endif


#endif