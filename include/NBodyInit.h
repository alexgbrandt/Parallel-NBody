
#ifndef _NBODY_INIT_
#define _NBODY_INIT_

#include <time.h>
#include <stdlib.h>

#include "NBodyHelpers.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Allocate the necessary arrays describing the n bodies.
 * All arrays returned are 3*n in size, except mass which is n.
 *
 * @param n, the number of bodies to allocate space for.
 * @param[out] r_p, return pointer for the position array.
 * @param[out] v_p, return pointer for the velocity array.
 * @param[out] a_p, return pointer for the acceleration array.
 * @param[out] m_p, return pointer for the mass array.
 * @param[out] work_p, return pointer for the work estiamtes array.
 * @return 0 iff the allocations was successful.
 */
int allocData(long n, double** r_p, double** v_p, double** a_p, double** m_p, double** work_p);


/**
 * Free the data created by allocData.
 *
 * @param r, the array of positions.
 * @param v, the array of velocities.
 * @param a, the array of acceleration.
 * @param m, the array of masses.
 */
static inline void freeData(double* r, double* v, double* a, double* m, double* work) {
    free(r);
    free(v);
    free(a);
    free(m);
    free(work);
}


/**
 * Initialize n masses so their masses are equal and sum is 1.0.
 *
 * @param n, the number of bodies.
 * @param m, an array of size n to store the masses.
 * @param M, the total mass of all bodies.
 *
 * @return 0 iff the initialization was successful.
 */
int _initMassEqual(long n, double* m, double M);


/**
 * Initialize n masses in solar mass units, based on the
 * mass generation function of Kroupa, Tout & Gilmore, 1993.
 * m(X) = 0.08 + g1*X^g2 + g3*X^g4 / (1-X)^0.58.
 * g1 = 0.19, g2 = 1.55, g3 = 0.05, g4=0.6.
 *
 * @param n, the number of bodies.
 * @param m, an array of size n to store the masses.
 *
 * @return 0 iff the initialization was successful.
 */
int _initMass(long n, double* m);


/**
 * Initialize n positions based on the Plummer model.
 * This method is based on Vol. 9 of Hut and Makino, 2005.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*n to store the positions.
 * @param Pr, Plummer radius
 *
 * @return 0 iff the initialization was successful.
 */
int _initPositionsPlummer(long n, double* r, double Pr);


/**
 * Initialize n positions uniformly throughout the unit sphere.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*n to store the positions.
 *
 * @return 0 iff the initialization was successful.
 */
int _initPositionsUniform(long n, double* r);


/**
 * Initialize n velocities based on the Plummer model
 * and Aaresh's escape velocity criteria.
 * This method is based on Vol. 9 of Hut and Makino, 2005.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*n which already stores the positions.
 * @param v, an array of size 3*n to store the velocities.
 *
 * @return 0 iff the initialization was successful.
 */
int _initVelocitiesPlummer(long n, const double* r, double* v);


/**
 * Initialize n velocities to zero.
 *
 * @param n, the number of bodies.
 * @param v, an array of size 3*n to store the velocities.
 *
 * @return 0 iff the initialization was successful.
 */
int _initVelocitiesZero(long n, double* v);


/**
 * Normalizes bodies into a reference frame where the center of mass
 * is at the origin with 0 velocity.
 *
 * @param n, the number of bodies
 * @param r, an array of 3*n doubles of positions.
 * @param v, an array of 3*n doubles of velocities.
 * @param m, an array of n doubles of masses.
 * @param M, the sum of masses.
 * @return 0 iff the adjustment was successful.
 */
int _centerOfMassAdjustment(long n, double* r, double* v, double* m, double M);


/**
 * Prepare the initial conditions of
 * position, velocity, acceleration, and mass
 * for n bodies.
 *
 * Normalizes bodies into a reference frame where the center of mass
 * is at the origin with 0 velocity.
 *
 * @param n, the number of bodies
 * @param seed, the random seed. If <= 0, use current time.
 * @param r, an array of 3*n doubles to store the positions.
 * @param v, an array of 3*n doubles to store the velocities.
 * @param a, an array of 3*n doubles to store the accelerations.
 * @param m, an array of n doubles to store the masses.
 *
 * @return 0 iff the initialization was successful.
 */
int initData(long n, time_t seed, double* r, double* v, double* a, double* m);


/**
 * Initialize the colors of the bodies for rendering.
 * These colors are based an interpolation of the apparent color
 * of main sequence stars of stellar type O through M.
 *
 * @param n, the number of colors to create.
 * @return an array of n RGBA colors as 4n floats.
 */
static inline float* createColors(long n) {
    //Polynomial fit coefs for color of main sequence stars, index equals monomial degree
    static float Rcoefs[5] = {1.0044296822813448, -0.34970836,  3.01217486, -7.03979628,  4.01061921};
    static float Gcoefs[5] = {0.7105430630750322, 0.35529946,  3.21427637, -8.17669968,  4.63094998};
    static float Bcoefs[6] = {0.47835976470247865, -0.88069571,  16.68565124, -44.11776039,  44.98967698, -16.16596219};
    float* color_data = (float*) malloc(sizeof(float)*4*n);
    float X[6];
    int k;
    for (long i = 0; i < n; ++i) {
        //X acts as an interpolation of spectral type between M and O as [0,1]
        X[1] = frand();
        for (k = 2; k < 6; ++k) {
            X[k] = X[k-1]*X[1];
        }
        color_data[4*i + 0] = Rcoefs[0] + X[1]*Rcoefs[1] + X[2]*Rcoefs[2] + X[3]*Rcoefs[3] + X[4]*Rcoefs[4];
        color_data[4*i + 1] = Gcoefs[0] + X[1]*Gcoefs[1] + X[2]*Gcoefs[2] + X[3]*Gcoefs[3] + X[4]*Gcoefs[4];
        color_data[4*i + 2] = Bcoefs[0] + X[1]*Bcoefs[1] + X[2]*Bcoefs[2] + X[3]*Bcoefs[3] + X[4]*Bcoefs[4] + X[5]*Bcoefs[5];
        color_data[4*i + 3] = 1.0;
    }
    return color_data;
}

#ifdef __cplusplus
}
#endif

#endif
