

#include "NBodyInit.h"
#include "NBodyForces.h"

#include "NBodyConfig.h"


/**
 * Allocate the necessary arrays describing the n bodies.
 * All arrays returned are 3*i in size, except mass which is n.
 *
 * @param n, the number of bodies to allocate space for.
 * @param[out] r_p, return pointer for the position array.
 * @param[out] v_p, return pointer for the velocity array.
 * @param[out] a_p, return pointer for the acceleration array.
 * @param[out] m_p, return pointer for the mass array.
 * @param[out] work_p, return pointer for the work estiamtes array.
 * @return 0 iff the allocations was successful.
 */
int allocData(long n, double** r_p, double** v_p, double** a_p, double** m_p, double** work_p) {
    //TODO consider using an alternating array here which
    //stores them all interleaved.
    if (r_p != NULL) {
        *r_p = (double*) malloc(sizeof(double)*n*3);
        if (*r_p == NULL) {
            return 1;
        }
    }
    if (v_p != NULL) {
        *v_p = (double*) malloc(sizeof(double)*n*3);
        if (*v_p == NULL) {
            return 1;
        }
    }
    if (a_p != NULL) {
        *a_p = (double*) malloc(sizeof(double)*n*3);
        if (*a_p == NULL) {
            return 1;
        }
    }
    if (m_p != NULL) {
        *m_p = (double*) malloc(sizeof(double)*n);
        if (*m_p == NULL) {
            return 1;
        }
    }
    if (work_p != NULL) {
        *work_p = (double*) malloc(sizeof(double)*n);
        if (*work_p == NULL) {
            return 1;
        }
    }

    return 0;
}


/**
 * Allocate an array to hold the work estimates for each particle.
 * These values are used for dynamic load balancing.
 *
 * @param n, the double of bodies to allocate for
 * @param[out] work_p, a pointer in which to store the allocated array.
 * @return 0 iff the allocations was successful.
 */
int allocWorkEstimates(long n, long** work_p) {
    if (work_p != NULL) {
        *work_p = (long*) malloc(sizeof(long)*n);
        if (*work_p == NULL) {
            return 1;
        }
    }

    return 0;
}


/**
 * Initialize n masses in solar mass units, based on the
 * mass generation function of Kroupa, Tout & Gilmore, 1993.
 * m(X) = 0.08 + g1*X^g2 + g3*X^g4 / (1-X)^0.58.
 * g1 = 0.19, g2 = 1.55, g3 = 0.05, g4=0.6.
 * Masses are scaled so that their sum is 1.0.
 *
 * @param n, the number of bodies.
 * @param m, an array of size n to store the masses.
 *
 * @return 0 iff the initialization was successful.
 */
int _initMass(long n, double* m) {
    const double g1 = 0.19;
    const double g2 = 1.55;
    const double g3 = 0.05;
    const double g4 = 0.6;

    double X, den;
    double M = 0.0;
    long i;
    for (i = 0; i < n; ++i) {
        X = frand();
        den = pow(1.0-X, 0.58);
        X = g1*pow(X, g2) + g3*pow(X, g4);
        X /= den;
        m[i] = X + 0.08;
        M += m[i];
    }

    M = 1.0/M;
    for (i = 0; i < n; ++i) {
        m[i] *= M;
    }

    return 0;
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
int _initMassEqual(long n, double* m, double M) {
    double mi = M / n;
    long i;

    for (i = 0; i < n; ++i) {
        m[i] = mi;
    }

    return 0;
}


/**
 * Initialize n positions based on the Plummer model.
 * This method is based on Vol. 9 of Hut and Makino, 2005.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*i to store the positions.
 * @param Pr, Plummer radius
 *
 * @return 0 iff the initialization was successful.
 */
int _initPositionsPlummer(long n, double* r, double Pr) {
    double R, X, Y;
    double n23 = -2.0 / 3.0;
    double Minv = 1.0; //assume we scale the mass so total mass is 1.0
    long i;
    Pr = Pr*Pr; //We only ever need it squared

    for (i = 0; i < n; ++i) {
        X = pow(frand()*Minv, n23);
        R = Pr / sqrt(X - 1);
        X = acos(1.0 - 2.0*frand()); //acos(-1..1) to get correct random dist.
        Y = 2.0*_PI*frand(); //phi

        r[3*i + 2] = R * cos(X); //z
        X = sin(X);
        r[3*i + 0] = R * X * cos(Y); //x
        r[3*i + 1] = R * X * sin(Y); //y
    }

    return 0;
}


/**
 * Initialize n positions uniformly throughout the unit sphere.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*n to store the positions.
 *
 * @return 0 iff the initialization was successful.
 */
int _initPositionsUniform(long n, double* r) {
    double R, X, Y;
    long i;
    for (i = 0; i < n; ++i) {
        R = frand();
        X = acos(1.0 - 2.0*frand());
        Y = frand()*2.0*_PI;
        r[3*i + 0] = R*sin(X)*cos(Y);
        r[3*i + 1] = R*sin(X)*sin(Y);
        r[3*i + 2] = R*cos(X);
    }

    return 0;
}


/**
 * Initialize n positions as two clusters separated by separation
 * in each dimension.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*n to store the positions.
 * @param separation, the distance separating the two clusters in each dimension.
 *
 * @return 0 iff the initialization was successful
 */
int _initTwoClusters(long n, double* r, double Pr, double* v, double separation) {
    long i = 0;
    long n1 = n/2;

    _initPositionsPlummer(n1, r, Pr);
    _initVelocitiesPlummer(n1, r, v);

    _initPositionsPlummer((n-n1), r + 3*n1, Pr);
    _initVelocitiesPlummer((n-n1), r+3*n1, v+3*n1);

    double dx = -separation, dy = -separation, dz = separation;
    for (i = 0; i < n1; ++i) {
        r[3*i] += dx;
        r[3*i + 1] += dy;
        r[3*i + 2] += dz;
    }

    dx = separation, dy = separation, dz = -separation;
    for (i = n1; i < n; ++i) {
        r[3*i] += dx;
        r[3*i + 1] += dy;
        r[3*i + 2] += dz;
    }

    return 0;
}


/**
 * Initialize n velocities based on the Plummer model
 * and Aaresh's escape velocity criteria.
 * This method is based on Vol. 9 of Hut and Makino, 2005.
 *
 * @param n, the number of bodies.
 * @param r, an array of size 3*i which already stores the positions.
 * @param v, an array of size 3*i to store the velocities.
 *
 * @return 0 iff the initialization was successful.
 */
int _initVelocitiesPlummer(long n, const double* r, double* v) {

    double R2, X, Y, vel;
    double sqrt2 = sqrt(2.0);
    long i;
    for (i = 0; i < n; ++i) {
        R2 = r[3*i]*r[3*i] + r[3*i + 1]*r[3*i + 1] + r[3*i + 2]*r[3*i + 2];

        do {
            X = frand();
            Y = 0.1*frand();
        } while (Y > X*X*pow(1.0-X*X, 3.5));
        vel = sqrt2 * X * pow(1 + R2, -0.25);

        X = acos(1.0 - 2.0*frand()); //acos(-1..1) to get correct random dist.
        Y = 2.0*_PI*frand(); //phi
        v[3*i + 2] = vel * cos(X); //z
        X = sin(X);
        v[3*i + 0] = vel * X * cos(Y); //x
        v[3*i + 1] = vel * X * sin(Y); //y
    }

    return 0;
}


/**
 * Initialize the n velocities randomly and uniformly in [-1,1] for each dimension.
 *
 * @param n, the number of velocities.
 * @param v, an array of 3*n doubles for n velocities.
 * @return 0 iff the initialization was successful.
 */
int _initVelocitiesUniform(long n, double* v) {
    long i;
    for (i = 0; i < n; ++i) {
        v[3*i] = (1.0 - 2.0*frand());
        v[3*i + 1] = (1.0 - 2.0*frand());
        v[3*i + 2] = (1.0 - 2.0*frand());
    }

    return 0;
}


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
int _centerOfMassAdjustment(long n, double* r, double* v, double* m, double M) {

    double rx = 0.0, ry = 0.0, rz = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double mi;
    long i;

    for (i = 0; i < n; ++i) {
        mi = m[i];
        rx += r[3*i]*mi;
        ry += r[3*i + 1]*mi;
        rz += r[3*i + 2]*mi;

        vx += v[3*i]*mi;
        vy += v[3*i + 1]*mi;
        vz += v[3*i + 2]*mi;
    }

    rx /= M;
    ry /= M;
    rz /= M;
    vx /= M;
    vy /= M;
    vz /= M;

    for (i = 0; i < n; ++i) {
        r[3*i] -= rx;
        r[3*i + 1] -= ry;
        r[3*i + 2] -= rz;
        v[3*i] -= vx;
        v[3*i + 1] -= vy;
        v[3*i + 2] -= vz;

    }

    return 0;
}


/**
 * Prepare the initial conditions of
 * position, velocity, acceleration, and mass
 * for n bodies.
 *
 * @param n, the number of bodies
 * @param seed, the random seed. If <= 0, use current time.
 * @param r, an array of 3*i doubles to store the positions.
 * @param v, an array of 3*i doubles to store the velocities.
 * @param a, an array of 3*i doubles to store the accelerations.
 * @param m, an array of n doubles to store the masses.
 *
 * @return 0 iff the initialization was successful.
 */
int initData(long n, time_t seed, double* r, double* v, double* a, double* m) {
    if (seed <= 0) {
        seed = time(NULL);
    }
    fprintf(stderr, "seed: %ld\n", seed);
    srand(seed);

    double M = 1.0;
    double Pr = 1.0;
    //scale factor so that total energy is -0.25 in Plummer distrib.
    double virialScale;
    virialScale = 16.0 / (3.0 * _PI); //for plummer model

    int error = _initMassEqual(n, m, M);
    if (error) {
        fprintf(stderr, "NBODY: Could not init masses.\n");
        exit(INIT_ERROR);
    }

    if (NBODY_INIT_TWO_CLUSTERS) {
        error = _initTwoClusters(n, r, Pr, v, 5.0);
    } else if (NBODY_INIT_PLUMMER) {
        // virialScale = 16.0 / (3.0 * _PI);
        error = _initPositionsPlummer(n, r, Pr);
        error = _initVelocitiesPlummer(n, r, v);
    } else if (NBODY_INIT_UNIFORM) {
        error = _initPositionsUniform(n, r);
        error = _initVelocitiesUniform(n, v);
    } else {
        //Plummer
        error = _initPositionsPlummer(n, r, Pr);
        error = _initVelocitiesPlummer(n, r, v);
    }
    if (error) {
        fprintf(stderr, "NBODY: Could not init nbody positions and velocities.\n");
        exit(INIT_ERROR);
    }

    error = _centerOfMassAdjustment(n, r, v, m, M);

    if (1) {
        //Aarseth, 2003, Algorithm 7.2.
        double Epot =  computeEpot_NB(n, m, r);
        double Ekin =  computeEkin_NB(n, m, v);
        double virialRatio = 0.5;
        double Qv = sqrt(virialRatio*fabs(Epot)/Ekin);
        scale3NArray_NB(n, v, Qv);
        double beta = (1 - virialRatio)*Epot/(Epot+Ekin);

        scale3NArray_NB(n, r, beta);
        scale3NArray_NB(n, v, 1.0/(sqrt(beta)));

        //After first scale Ekin is -0.5Epot but E0 != -0.25.
        //So just scale up or down as needed.
        Epot = computeEpot_NB(n, m, r);
        beta = Epot / -0.5;
        scale3NArray_NB(n, r, beta);
        scale3NArray_NB(n, v, 1.0/sqrt(beta));
    } else {
        scale3NArray_NB(n, r, 1.0/virialScale);
        scale3NArray_NB(n, v, sqrt(virialScale));

    }

    if (error) {
        fprintf(stderr, "NBODY: Could not scale to standard units.\n");
        exit(INIT_ERROR);
    }

    computeForces(n, m, r, a);

    long i;
    for (i = 0; i < n; ++i) {
        a[3*i + 0] = 0.0;
        a[3*i + 1] = 0.0;
        a[3*i + 2] = 0.0;
    }

    return (error == 0);
}


/**
 * Chapter 5: Art of computational science vol_1 v1_web.
 */
int initData3BodyChaotic(double* r, double* v, double* a, double* m) {
    m[0] = 1.0;
    m[1] = 1.0;
    m[2] = 1.0;
    long i;
    for (i = 0; i < 3; ++i) {
        double phi = i*2*_PI / 3.0;
        r[3*i + 0] = cos(phi);
        r[3*i + 1] = sin(phi);
        r[3*i + 2] = 0.0;
    }

    computeForces(3, m, r, a);

    double v_abs = sqrt(-1.0*a[0]);
    for (i = 0; i < 3; ++i) {
        double phi = i*2*_PI / 3.0;
        v[3*i + 0] = -1.0 * v_abs * sin(phi);
        v[3*i + 1] = v_abs * cos(phi);
        v[3*i + 2] = 0.0;
    }
    v[0] += 0.0001;

    return 0;
}


/**
 * 3Body figure eight stable orbit init.
 */
int initData3BodyFigureEight(double* r, double* v, double* a, double* m) {

    m[0] = 1.0;
    m[1] = 1.0;
    m[2] = 1.0;

    r[3*0 + 0] = 0.9700436;
    r[3*0 + 1] = -0.24308753;
    r[3*0 + 2] = 0.0;
    v[3*0 + 0] = 0.466203685;
    v[3*0 + 1] = 0.43236573;
    v[3*0 + 2] = 0.0;

    r[3*1 + 0] = -r[3*0 + 0];
    r[3*1 + 1] = -r[3*0 + 1];
    r[3*1 + 2] = -r[3*0 + 2];
    v[3*1 + 0] = v[3*0 + 0];
    v[3*1 + 1] = v[3*0 + 1];
    v[3*1 + 2] = v[3*0 + 2];

    r[3*2 + 0] = 0.0;
    r[3*2 + 1] = 0.0;
    r[3*2 + 2] = 0.0;

    v[3*2 + 0] = -2.0 * v[3*0 + 0];
    v[3*2 + 1] = -2.0 * v[3*0 + 1];
    v[3*2 + 2] = -2.0 * v[3*0 + 2];

    v[0] += 0.001;

    computeForces(3, m, r, a);

    return 0;
}
