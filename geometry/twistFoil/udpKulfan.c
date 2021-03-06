/*
 ************************************************************************
 *                                                                      *
 * udpKulfan -- udp file to generate Kulfan airfoils                    *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 *            Based upon prescription by B.M. Kulfan and J.E. Bussoletti*
 *              "Fundamental Parametric Geometry Representations for    *
 *               Aircraft Component Shapes", AIAA-2006-6948             *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2022  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#include "egads_dot.h"

#define NUMUDPARGS 5
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define CLASS(     IUDP,I)  ((double *) (udps[IUDP].arg[0].val))[I]
#define CLASS_DOT( IUDP,I)  ((double *) (udps[IUDP].arg[0].dot))[I]
#define ZTAIL(     IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define ZTAIL_DOT( IUDP,I)  ((double *) (udps[IUDP].arg[1].dot))[I]
#define AUPPER(    IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]
#define AUPPER_DOT(IUDP,I)  ((double *) (udps[IUDP].arg[2].dot))[I]
#define ALOWER(    IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define ALOWER_DOT(IUDP,I)  ((double *) (udps[IUDP].arg[3].dot))[I]
#define NUMPTS    (IUDP  )  ((int    *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"class",     "ztail",     "aupper",    "alower",    "numpts", };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, ATTRINT,  };
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,           0,           101,      };
static double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,          0.,          0.,       };

typedef struct {
  int    sharpte;
  double class_dot[2];
  double ztail_dot[2];
  double *aupper_dot;
  double *alower_dot;
} udpDotCache_T;

#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data)
{
  udpDotCache_T *cache = (udpDotCache_T*)data;
  if (cache != NULL) {
      EG_free(cache->aupper_dot);
      EG_free(cache->alower_dot);
      EG_free(cache);
  }
  return EGADS_SUCCESS;
}

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
 udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

/* set KNOTS to 0 for arc-lenght knots, and -1 for equally spaced knots
 */
#define           KNOTS           0

#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           DXYTOL          1.0e-6
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

#ifdef GRAFIC
    #include "grafic.h"
#endif

/*
 ************************************************************************
 *                                                                      *
 *   factorial - factorial of a number (with cache)                     *
 *                                                                      *
 ************************************************************************
 */

static long
factorial(int n)
{
    long result;
    static long table[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                           362880, 3628800, 39916800, 479001600};

    /* bad input */
    if        (n < 0) {
        result = 0;

    /* use (cached) result from table */
    } else if (n < 13) {
        result = table[n];

    /* compute the result if not in the table */
    } else {
        result = n;
        while (--n > 1) {
            result *= n;
        }
    }

    return result;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     i, r, n, ipnt, sense[3], nedge, *header=NULL, sizes[2], oclass, mtype;
    double  zeta, shape, s, K, pow1, pow2, *pnts=NULL, *rdata=NULL;
    double  data[18], tdata[2], tle;
    ego     eref, enodes[4], eedges[3], ecurve, eline, eloop, eplane, eface;
    udpDotCache_T *cache=NULL;
#ifdef GRAFIC
    float   *xplot=NULL, *yplot=NULL;
#endif

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if        (udps[0].arg[0].size < 2) {
        printf(" udpExecute: class should contain 2 values (nose,tail)\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size < 2) {
        printf(" udpExecute: ztail should contain 2 values (upper,lower)\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size < 2) {
        printf(" udpExecute: aupper should contain at least 1 value\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size < 2) {
        printf(" udpExecute: aupper should contain at least 1 value\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size != 1 || NUMPTS(0) < 11) {
        printf(" udpExecute: numpts should contain one positive number greater than 11\n");
        status = EGADS_RANGERR;
        goto cleanup;

    }

#ifdef DEBUG
    {
        int i;

        printf("class( 0)[0]   = %f  (nose)\n", CLASS(0,0));
        printf("class( 0)[1]   = %f  (tail)\n", CLASS(0,1));
        printf("ztail( 0)[0]   = %f  (upper)\n", ZTAIL(0,0));
        printf("ztail( 0)[1]   = %f  (lower)\n", ZTAIL(0,1));
        for (i = 0; i < udps[0].arg[2].size; i++) {
            printf("aupper(0)[%d]  = %f\n", i, AUPPER(0,i));
        }
        for (i = 0; i < udps[0].arg[3].size; i++) {
            printf("alower(0)[%d]  = %f\n", i, ALOWER(0,i));
        }
        printf("numpts(0)      = %d\n", NUMPTS( 0));
    }
#endif

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

    cache = (udpDotCache_T*)EG_alloc(sizeof(udpDotCache_T));
    udps[numUdp].data = (void*)cache;

    cache->sharpte      = 1;
    cache->class_dot[0] = 0;
    cache->class_dot[1] = 0;
    cache->ztail_dot[0] = 0;
    cache->ztail_dot[1] = 0;
    cache->aupper_dot = (double*)EG_alloc(udps[numUdp].arg[2].size*sizeof(double));
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
      cache->aupper_dot[i] = 0.0;
    }
    cache->alower_dot = (double*)EG_alloc(udps[numUdp].arg[3].size*sizeof(double));
    for (i = 0; i < udps[numUdp].arg[3].size; i++) {
      cache->aupper_dot[i] = 0.0;
    }

    /* mallocs required by Windows compiler */
    MALLOC(pnts, double, (3*NUMPTS(numUdp)));

    /* points around airfoil ( upper and lower) */
    for (ipnt = 0; ipnt < NUMPTS(numUdp); ipnt++) {
        zeta = TWOPI * ipnt / (NUMPTS(numUdp)-1);
        s    = (1 + cos(zeta)) / 2;

        /* upper surface */
        if ( ipnt < (NUMPTS(numUdp)-1)/2){
            n = udps[numUdp].arg[2].size - 1;

            shape = 0;
            for (r = 0; r <= n; r++) {
                pow1 = pow(1-s, n-r);
                pow2 = pow(  s,   r);
                K    = factorial(n) / factorial(r) / factorial(n-r);
                shape += AUPPER(numUdp,r) * K * pow1 * pow2;
            }

            pow1 = pow(  s, CLASS(numUdp,0));
            pow2 = pow(1-s, CLASS(numUdp,1));
            pnts[3*ipnt  ] = s;
            pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(numUdp,0) * s;
            pnts[3*ipnt+2] = 0;

        /* leading edge */
        } else if (ipnt == (NUMPTS(numUdp)-1)/2) {
            pnts[3*ipnt  ] = 0;
            pnts[3*ipnt+1] = 0;
            pnts[3*ipnt+2] = 0;

        /* lower surface */
        } else if (ipnt > (NUMPTS(numUdp)-1)/2) {
            n = udps[numUdp].arg[3].size - 1;

            shape = 0;
            for (r = 0; r <= n; r++) {
                pow1 = pow(1-s, n-r);
                pow2 = pow(  s,   r);
                K = factorial(n) / factorial(r) / factorial(n-r);
                shape += ALOWER(numUdp,r) * K * pow1 * pow2;
            }

            pow1 = pow(  s, CLASS(numUdp,0));
            pow2 = pow(1-s, CLASS(numUdp,1));
            pnts[3*ipnt  ] = s;
            pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(numUdp,1) * s;
            pnts[3*ipnt+2] = 0;
        }
    }

    /* create spline curve
     *
     * equally spaced (sizes[1] == -1)
     * arc-length based knots (sizes[1] == 0)
     */
    sizes[0] = NUMPTS(numUdp);
    sizes[1] = KNOTS;
    status = EG_approximate(context, 0, DXYTOL, sizes, pnts, &ecurve);
    CHECK_STATUS(EG_approximate);

#ifdef GRAFIC
    #define NUMEVAL 5000

    if (1) {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    nline=0, npnt=0, ilin[NUMEVAL+3], isym[NUMEVAL+3], nper[NUMEVAL+3], i;
        double tt, evaldata[18];

        MALLOC(xplot, float, NUMPTS(numUdp)+3*NUMEVAL);
        MALLOC(yplot, float, NUMPTS(numUdp)+3*NUMEVAL);

        for (i = 0; i < NUMPTS(numUdp); i++) {
            xplot[npnt] = pnts[3*i  ];
            yplot[npnt] = pnts[3*i+1];
            npnt++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] =  NUMPTS(numUdp);
        nline++;

        for (i = 0; i < NUMEVAL; i++) {
            tt = (double)(i) / (double)(NUMEVAL-1);
            (void) EG_evaluate(ecurve, &tt, evaldata);

            xplot[npnt] = evaldata[0];
            yplot[npnt] = evaldata[1];
            npnt++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] =  NUMEVAL;
        nline++;

        for (i = 0; i < NUMEVAL; i++) {
            tt = (double)(i) / (double)(NUMEVAL+1);
            (void) EG_evaluate(ecurve, &tt, evaldata);

            xplot[npnt] = evaldata[0];
            yplot[npnt] = evaldata[1];
            npnt++;

            xplot[npnt] = evaldata[0] + evaldata[4] / 100;
            yplot[npnt] = evaldata[1] - evaldata[3] / 100;
            npnt++;

            ilin[nline] = GR_SOLID;
            isym[nline] = 0;
            nper[nline] = 2;
            nline++;
        }

        grinit_(&io_kbd, &io_scr, "udpKulfan", STRLEN("udpKulfan"));
        grline_(ilin, isym, &nline,                 "~x~y~O=points, line=fit",
                &indgr, xplot, yplot, nper, STRLEN( "~x~y~O=points, line=fit"));
    }
#endif

    /* create Node at upper trailing edge */
    ipnt = 0;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    /* Node at leading edge as a function of the spline */
    status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &header, &rdata);
    CHECK_STATUS(EG_getGeometry);

    ipnt = (NUMPTS(numUdp) - 1) / 2 + 3; /* index, with knot offset of 3 (cubic)*/
    tle  = rdata[ipnt];           /* t-value (should be very close to (0,0,0) */

    status = EG_evaluate(ecurve, &tle, data);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* make Edge for upper surface */
    tdata[0] = 0;   /* t-value at lower TE */
    tdata[1] = tle;

    /* construct the upper Edge */
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    /* create line segment at trailing edge */
    ipnt = NUMPTS(numUdp) - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    data[3] = pnts[0] - data[0];
    data[4] = pnts[1] - data[1];
    data[5] = pnts[2] - data[2];

    if (fabs(data[3]) > EPS06 || fabs(data[4]) > EPS06 || fabs(data[5]) > EPS06) {
        nedge = 3;
        cache->sharpte = 0;

        /* create Node at lower trailing edge */
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enodes[2]));
        CHECK_STATUS(EG_makeTopology);

        /* make Edge for lower surface */
        tdata[0] = tdata[1]; /* t-value at leading edge */
        tdata[1] = 1;        /* t value at upper TE */

        /* construct the lower Edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &enodes[1], NULL, &eedges[1]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
        CHECK_STATUS(EG_makeGeometry);

        /* reuse the trailing edge node */
        enodes[3] = enodes[0];

        /* make Edge for this line */
        tdata[0] = 0;
        tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

        status = EG_makeTopology(context, eline, EDGE, TWONODE,
                                 tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
        CHECK_STATUS(EG_makeTopology);
    } else {
        nedge = 2;
        cache->sharpte = 1;

        /* reuse the trailing edge node */
        enodes[2] = enodes[0];

        /* make Edge for lower surface */
        tdata[0] = tdata[1]; /* t-value at leading edge */
        tdata[1] = 1;        /* t value at upper TE */

        /* construct the lower Edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);
    }

    /* create loop of either two or three Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, nedge, eedges, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    /* create a plane for the loop */
    data[0] = 0.;
    data[1] = 0.;
    data[2] = 0.;
    data[3] = 1.; data[4] = 0.; data[5] = 0.;
    data[6] = 0.; data[7] = 1.; data[8] = 0.;

    status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, data, &eplane);
    CHECK_STATUS(EG_makeGeometry);

    /* create the face from the plane and the loop */
    status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                             NULL, 1, &eloop, sense, &eface);
    CHECK_STATUS(EG_makeTopology);

    /* create the face body */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* set the output value(s) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    EG_free(header);
    EG_free(rdata);
    EG_free(pnts);

#ifdef GRAFIC
    FREE(xplot);
    FREE(yplot);
#endif

    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }

    return status;
}

/*
 *****************************************************************************
 *                                                                           *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                           *
 *****************************************************************************
 */

int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp;

    int    ipnt, nedge, i, r, n, idotChange, sizes[2];
    int    oclass, mtype, nchild, *senses, stride;
    double shape, shape_dot, s, K, pow1, pow2, zeta;
    double *pnts=NULL, *pnts_dot=NULL, point[18], point_dot[18];
    double  data[18], data_dot[18], tdata[2], tdata_dot[2];
    double  tle, tle_dot, *rdata=NULL, *rdata_dot=NULL;;
    ego    eref, *echildren, eface, eplane, eloop, enodes[3], *eedges, ecurve, eline, eent;
    udpDotCache_T *cache;

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
#endif

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    /* check if velocities have changed */
    cache = (udpDotCache_T*)udps[iudp].data;
    idotChange = 0;
    for (i = 0; i < 2; i++) {
        if ( cache->class_dot[i] != CLASS_DOT(iudp, i) ) {
            idotChange = 1;
        }
        if ( cache->ztail_dot[i] != ZTAIL_DOT(iudp, i) ) {
            idotChange = 1;
        }
    }
    for (i = 0; i < udps[iudp].arg[2].size && idotChange == 0; i++)
        if ( cache->aupper_dot[i] != AUPPER_DOT(iudp,i) ) {
            idotChange = 1;
        }
    for (i = 0; i < udps[iudp].arg[3].size && idotChange == 0; i++)
        if ( cache->alower_dot[i] != ALOWER_DOT(iudp,i) ) {
            idotChange = 1;
        }

    /* populate the egos with sensitivities if velocities changed or do not exist */
    if (idotChange == 1 ||
        EG_hasGeometry_dot(ebody) != EGADS_SUCCESS) {

        for (i = 0; i < 2; i++) {
            cache->class_dot[i] = CLASS_DOT(iudp, i);
            cache->ztail_dot[i] = ZTAIL_DOT(iudp, i);
        }
        for (i = 0; i < udps[iudp].arg[2].size; i++) {
            cache->aupper_dot[i] = AUPPER_DOT(iudp,i);
        }
        for (i = 0; i < udps[iudp].arg[3].size; i++) {
            cache->alower_dot[i] = ALOWER_DOT(iudp,i);
        }

        /* get the face from the FACEBODY */
        status = EG_getTopology(ebody, &eref, &oclass, &mtype, data, &nchild, &echildren,
                                &senses);
        CHECK_STATUS(EG_getTopology);
        eface = echildren[0];

        /* get the plane and the loop */
        status = EG_getTopology(eface, &eplane, &oclass, &mtype, data, &nchild, &echildren,
                                &senses);
        CHECK_STATUS(EG_getTopology);
        eloop = echildren[0];

        /* get the edges from the loop */
        status = EG_getTopology(eloop, &eref, &oclass, &mtype, data, &nedge, &eedges,
                                &senses);
        CHECK_STATUS(EG_getTopology);

        /* get the nodes and the curve from the first edge */
        status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, data, &nchild, &echildren,
                                &senses);
        CHECK_STATUS(EG_getTopology);
        enodes[0] = echildren[0]; /* upper trailing edge */
        enodes[1] = echildren[1]; /* leading edge        */

        /* create the kulfan spline points with sensitivities */

        /* mallocs required by Windows compiler */
        MALLOC(pnts    , double, (3*NUMPTS(iudp)));
        MALLOC(pnts_dot, double, (3*NUMPTS(iudp)));

        /* points around airfoil ( upper and lower) */
        for (ipnt = 0; ipnt < NUMPTS(iudp); ipnt++) {
            zeta = TWOPI * ipnt / (NUMPTS(iudp)-1);
            s    = (1 + cos(zeta)) / 2;

            /* upper surface */
            if ( ipnt < (NUMPTS(iudp)-1)/2){
                n = udps[iudp].arg[2].size - 1;

                shape     = 0;
                shape_dot = 0;
                for (r = 0; r <= n; r++) {
                    pow1 = pow(1-s, n-r);
                    pow2 = pow(  s,   r);
                    K    = factorial(n) / factorial(r) / factorial(n-r);
                    shape     += AUPPER(    iudp,r) * K * pow1 * pow2;
                    shape_dot += AUPPER_DOT(iudp,r) * K * pow1 * pow2;
                }

                /* points for the spline fit */
                pow1 = pow(  s, CLASS(iudp,0));
                pow2 = pow(1-s, CLASS(iudp,1));
                pnts[3*ipnt  ] = s;
                pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(iudp,0) * s;
                pnts[3*ipnt+2] = 0;

                /* velocity of the points for the spline fit */
                pnts_dot[3*ipnt  ] = 0;
                pnts_dot[3*ipnt+1] = pow1 * pow2 * shape_dot + ZTAIL_DOT(iudp,0) * s;
                if (s   > 0) pnts_dot[3*ipnt+1] += pow1 * pow2 * shape * log(  s)
                                                   * CLASS_DOT(iudp,0);
                if (1-s > 0) pnts_dot[3*ipnt+1] += pow1 * pow2 * shape * log(1-s)
                                                   * CLASS_DOT(iudp,1);
                pnts_dot[3*ipnt+2] = 0;

            /* leading edge */
            } else if (ipnt == (NUMPTS(iudp)-1)/2) {
                pnts[3*ipnt  ] = 0;
                pnts[3*ipnt+1] = 0;
                pnts[3*ipnt+2] = 0;

                pnts_dot[3*ipnt  ] = 0;
                pnts_dot[3*ipnt+1] = 0;
                pnts_dot[3*ipnt+2] = 0;

            /* lower surface */
            } else if (ipnt > (NUMPTS(iudp)-1)/2) {
                n = udps[iudp].arg[3].size - 1;

                shape     = 0;
                shape_dot = 0;
                for (r = 0; r <= n; r++) {
                    pow1 = pow(1-s, n-r);
                    pow2 = pow(  s,   r);
                    K    = factorial(n) / factorial(r) / factorial(n-r);
                    shape     += ALOWER(    iudp,r) * K * pow1 * pow2;
                    shape_dot += ALOWER_DOT(iudp,r) * K * pow1 * pow2;
                }

                /* points for the spline fit */
                pow1 = pow(  s, CLASS(iudp,0));
                pow2 = pow(1-s, CLASS(iudp,1));
                pnts[3*ipnt  ] = s;
                pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(iudp,1) * s;
                pnts[3*ipnt+2] = 0;

                /* velocity of the points for the spline fit */
                pnts_dot[3*ipnt  ] = 0;
                pnts_dot[3*ipnt+1] = pow1 * pow2 * shape_dot + ZTAIL_DOT(iudp,1) * s;
                if (s   > 0) pnts_dot[3*ipnt+1] += pow1 * pow2 * shape * log(  s)
                                                   * CLASS_DOT(iudp,0);
                if (1-s > 0) pnts_dot[3*ipnt+1] += pow1 * pow2 * shape * log(1-s)
                                                   * CLASS_DOT(iudp,1);
                pnts_dot[3*ipnt+2] = 0;
            }
        }

        /* populate spline curve sensitivities */
        sizes[0] = NUMPTS(iudp);
        sizes[1] = KNOTS;
        status = EG_approximate_dot(ecurve, 0, DXYTOL, sizes, pnts, pnts_dot);
        CHECK_STATUS(EG_approximate_dot);

        /* set the sensitivity of the Node at trailing edge */
        ipnt = 0;
        status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, &(pnts[3*ipnt]), &(pnts_dot[3*ipnt]));
        CHECK_STATUS(EG_setGeometry_dot);


        /* set the sensitivity of the Node at leading edge */
        status = EG_getGeometry_dot(ecurve, &rdata, &rdata_dot);
        CHECK_STATUS(EG_getGeometry_dot);

        ipnt = (NUMPTS(iudp) - 1) / 2 + 3; /* index, with knot offset of 3 (cubic)*/
        tle     = rdata[ipnt];        /* t-value (should be very close to (0,0,0) */
        tle_dot = rdata_dot[ipnt];    /* t-value sensitivity */

        status = EG_evaluate_dot(ecurve, &tle, &tle_dot, data, data_dot);
        CHECK_STATUS(EG_evaluate_dot);
        status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
        CHECK_STATUS(EG_setGeometry_dot);


        /* set Edge t-range sensitivity for upper surface */
        tdata[0]     = 0;
        tdata[1]     = tle;
        tdata_dot[0] = 0;
        tdata_dot[1] = tle_dot;

        status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
        CHECK_STATUS(EG_setRange_dot);

        /* set Edge t-range sensitivity for lower surface */
        tdata[0]     = tdata[1];
        tdata[1]     = 1;
        tdata_dot[0] = tdata_dot[1];
        tdata_dot[1] = 0;

        status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
        CHECK_STATUS(EG_setRange_dot);


        if (cache->sharpte == 0) {
            /* get trailing edge line and the lower trailing edge node from the 3rd edge */
            status = EG_getTopology(eedges[2], &eline, &oclass, &mtype, data, &nchild, &echildren,
                                    &senses);
            CHECK_STATUS(EG_getTopology);
            enodes[2] = echildren[0]; /* lower trailing edge */

            /* set the sensitivity of the Node at lower trailing edge */
            ipnt = NUMPTS(iudp) - 1;
            status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, &(pnts[3*ipnt]), &(pnts_dot[3*ipnt]));
            CHECK_STATUS(EG_setGeometry_dot);

            /* set the sensitivity of the line segment at trailing edge */
            ipnt = NUMPTS(iudp) - 1;
            data[0] = pnts[3*ipnt  ];
            data[1] = pnts[3*ipnt+1];
            data[2] = pnts[3*ipnt+2];
            data[3] = pnts[0] - data[0];
            data[4] = pnts[1] - data[1];
            data[5] = pnts[2] - data[2];

            data_dot[0] = pnts_dot[3*ipnt  ];
            data_dot[1] = pnts_dot[3*ipnt+1];
            data_dot[2] = pnts_dot[3*ipnt+2];
            data_dot[3] = pnts_dot[0] - data_dot[0];
            data_dot[4] = pnts_dot[1] - data_dot[1];
            data_dot[5] = pnts_dot[2] - data_dot[2];

            status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
            CHECK_STATUS(EG_setGeometry_dot);

            /* set Edge t-range sensitivity */
            tdata[0] = 0;
            tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

            tdata_dot[0] = 0;
            tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

            status = EG_setRange_dot(eedges[2], EDGE, tdata, tdata_dot);
            CHECK_STATUS(EG_setRange_dot);
        }

        /* plane data */
        data[0] = 0.;
        data[1] = 0.;
        data[2] = 0.;
        data[3] = 1.; data[4] = 0.; data[5] = 0.;
        data[6] = 0.; data[7] = 1.; data[8] = 0.;

        /* set the sensitivity of the plane */
        data_dot[0] = 0.;
        data_dot[1] = 0.;
        data_dot[2] = 0.;
        data_dot[3] = 0.; data_dot[4] = 0.; data_dot[5] = 0.;
        data_dot[6] = 0.; data_dot[7] = 0.; data_dot[8] = 0.;

        status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, data, data_dot);
        CHECK_STATUS(EG_setGeometry_dot);

    } /* if (udps[iudp].data == NULL) */


    /* find the ego entity */
    if (entType == OCSM_NODE) {
        status = EG_getBodyTopos(ebody, NULL, NODE, &nchild, &echildren);
        CHECK_STATUS(EG_getBodyTopos);

        stride = 0;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else if (entType == OCSM_EDGE) {
        status = EG_getBodyTopos(ebody, NULL, EDGE, &nchild, &echildren);
        CHECK_STATUS(EG_getBodyTopos);

        stride = 1;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else if (entType == OCSM_FACE) {
        status = EG_getBodyTopos(ebody, NULL, FACE, &nchild, &echildren);
        CHECK_STATUS(EG_getBodyTopos);

        stride = 2;
        eent = echildren[entIndex-1];

        EG_free(echildren); echildren = NULL;
    } else {
        printf("udpSensitivity: bad entType=%d\n", entType);
        status = EGADS_GEOMERR;
        goto cleanup;
    }

    /* get the velocities from the entity */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        status = EG_evaluate_dot(eent, &(uvs[stride*ipnt]), NULL, point, point_dot);
        CHECK_STATUS(EG_evaluate_dot);

        /* return the point velocity */
        vels[3*ipnt  ] = point_dot[0];
        vels[3*ipnt+1] = point_dot[1];
        vels[3*ipnt+2] = point_dot[2];
    }

    status = EGADS_SUCCESS;

cleanup:
    EG_free(pnts);
    EG_free(pnts_dot);
    EG_free(rdata);
    EG_free(rdata_dot);

    return status;
}
