/*

 $Id: qmultiply.cc,v 1.15 2009/05/08 23:02:16 rhuey Exp $
 $Id: qmultiply.cc, autodock4.lga.MPI v0.0 2010/10/10 22:55:01 collignon Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "qmultiply.h"

extern  FILE    *logFile;


void qmultiply( Quat *q, //result
                register const Quat *ql,   //left
                register const Quat *qr )  //right

/******************************************************************************/
/*      Name: qmultiply                                                       */
/*  Function: Quaternion Multiplication (Accelerated)                         */
/*            [q]  =  [ql] [qr]                                               */
/*            [s1,v1][s2,v2] = [(s1*s2 - v1.v2), (s1*v2 + s2*v1 + v1^v2)]     */
/*                ~~     ~~              ~~ ~~       ~~      ~~   ~~ ~~       */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, The Scripps Research Institute.              */
/*            David Goodsell, TSRI                                            */
/*      Date: 12/03/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: ql = rotation to be applied to quaternion in qr                 */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/15/92 GMM     Translated into C                                         */
/* 12/03/92 GMM     Changed '/2.' to '*0.5'; introduced 'hqwl' and 'hqwr'.    */
/* 12/03/92 GMM     Replaced rqtot by inv_qag; was '/rqtot', now '*inv_qag'   */
/******************************************************************************/
{ 
    register double x,y,z,w;

    x = (double) (ql->w*qr->x + ql->x*qr->w + ql->y*qr->z - ql->z*qr->y);
    y = (double) (ql->w*qr->y + ql->y*qr->w + ql->z*qr->x - ql->x*qr->z);
    z = (double) (ql->w*qr->z + ql->z*qr->w + ql->x*qr->y - ql->y*qr->x);
    w = (double) (ql->w*qr->w - ql->x*qr->x - ql->y*qr->y - ql->z*qr->z);

    q->x = x;
    q->y = y;
    q->z = z;
    q->w = w;
}

void qconjmultiply( Quat *q,
                    register const Quat *ql,
                    register const Quat *qr )
//     __     
// q = ql . qr
{
    Quat conj_ql = conjugate( *ql );
    qmultiply( q, &conj_ql, qr );
}

void mkUnitQuat( Quat *q )
    // essentially, convertRotToQuat( Quat q )
{	
    double inv_nmag, hqang, s;
	     
    inv_nmag = 1. / hypotenuse( q->nx, q->ny, q->nz );
    q->nx *= inv_nmag;       /* Normalize q */
    q->ny *= inv_nmag;       /* Normalize q */
    q->nz *= inv_nmag;       /* Normalize q */
      
    hqang = 0.5 * q->ang;
    s     = sin( hqang );
    
    q->w  = cos( hqang );
    q->x  = s * q->nx;
    q->y  = s * q->ny;
    q->z  = s * q->nz;
    
    /* q->qmag = hypotenuse4( q->x,  q->y,  q->z,  q->w  ); */
} // mkUnitQuat( Quat *q )

void printQuat_q( FILE *fp, Quat q )
{
    (void) fprintf( fp, "Quat(x,y,z,w)=        %5.2f %5.2f %5.2f %5.2f\n", q.x, q.y, q.z, q.w);
    (void) fprintf( fp, "Mag(Quat(x,y,z,w))=   %5.2f\n", sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w) );
} // printQuat_q( Quat q )

void printQuat_r( FILE *fp, Quat q )
{
    q = convertQuatToRot( q );
    (void) fprintf( fp, "Axis(nx,ny,nz),Angle= %5.2f %5.2f %5.2f  %5.2f\n", q.nx, q.ny, q.nz, q.ang);
    (void) fprintf( fp, "Mag(Axis(nx,ny,nz))=  %5.2f\n", sqrt(q.nx*q.nx + q.ny*q.ny + q.nz*q.nz) );
} // printQuat_r( Quat q )

void printQuat( FILE *fp, Quat q )
{
    printQuat_q( fp, q );
    printQuat_r( fp, q );
} // printQuat( Quat q )

void debugQuat( FILE *fp, Quat q, unsigned int linenumber, char *message )
{
    pr( fp, "DEBUG_QUAT: %s   (line %u)\n", message,  linenumber );
    printQuat( fp, q );
}

Quat normQuat( Quat q )
    // Normalise the 4D quaternion, x,y,z,w
{
    register double mag4 = hypotenuse4( q.x, q.y, q.z, q.w );
    if (mag4 > APPROX_ZERO) {
        register double inv_mag4 = 1. / mag4;
        q.x *= inv_mag4;
        q.y *= inv_mag4;
        q.z *= inv_mag4;
        q.w *= inv_mag4;
    }
    return q;
}

/*
#define ONE_MINUS_EPSILON 0.999
#define ONE_PLUS_EPSILON 1.001

void assertQuatOK( const Quat q )
{
    register double mag4 = hypotenuse4( q.x, q.y, q.z, q.w );
    assert((mag4 > ONE_MINUS_EPSILON) && (mag4 < ONE_PLUS_EPSILON));
}
*/
/* this is in another header
#define assertQuatOK( q ) {register double aQOK_mag4 = hypotenuse4( (q).x, (q).y, (q).z, (q).w ); assert((aQOK_mag4 > ONE_MINUS_EPSILON) && (aQOK_mag4 < ONE_PLUS_EPSILON)); }
*/

Quat normRot( Quat q )
    // Normalise the 3D rotation axis or vector nx,ny,nz
{
    double mag3 = hypotenuse( q.nx, q.ny, q.nz );
    if (mag3 > APPROX_ZERO) {
        double inv_mag3 = 1. / mag3;
        q.nx *= inv_mag3;
        q.ny *= inv_mag3;
        q.nz *= inv_mag3;
    }
    return q;
}

Real quatDifferenceToAngle( const Quat ql, const Quat qr )
{
    Quat qdiff, rot;
    qconjmultiply(&qdiff, &ql, &qr);
    rot = convertQuatToRot( qdiff );
    return rot.ang;
}

Real quatDifferenceToAngleDeg( const Quat ql, const Quat qr )
{
    return (180./PI)* quatDifferenceToAngle( ql, qr );
}


Quat convertQuatToRot( Quat q )
    // Convert the quaternion components (x,y,z,w) of the quaternion q,
    // to the corresponding rotation-about-axis components (nx,ny,nz,ang)
{
    // TODO handle big W!  Singularities...
    Quat retval;
#ifdef SUPER_DEBUG_MUTATION // mp
    fprintf( logFile, "convertQuatToRot:  q.w = %.3f\n", q.w );
#endif
    assert( fabs( q.w ) <= 1.001 );
    if ( q.w > 1. ) q.w = 1.;
    if ( q.w < -1. ) q.w = -1.;

    register double angle = 2. * acos( q.w );
    register double inv_sin_half_angle = 1.;
    if ( q.w == 1. ) {
        retval.nx = 1.;
        retval.ny = 0.;
        retval.nz = 0.;
    } else {
        inv_sin_half_angle = 1. / sin( angle / 2. );

        retval.nx = q.x * inv_sin_half_angle;
        retval.ny = q.y * inv_sin_half_angle;
        retval.nz = q.z * inv_sin_half_angle;

        retval = normRot( retval );
    }
    angle = WrpModRad( angle );  // by convention, angles should be in the range -PI to +PI.
    retval.ang = angle;

    // Copy the existing x,y,z,w components
    retval.x = q.x;
    retval.y = q.y;
    retval.z = q.z;
    retval.w = q.w;

    return retval;
} // convertQuatToRot( Quat q )

Quat convertRotToQuat( Quat q )
    // Normalize the rotation-about-axis vector 
    // and convert the rotation-about-axis components (nx,ny,nz,ang)
    // to the corresponding quaternion components (x,y,z,w)
{	
    double hqang, s;
    Quat retval;

    retval.nx = q.nx;
    retval.ny = q.ny;
    retval.nz = q.nz;
    retval = normRot( retval );

    retval.ang = q.ang;
      
    hqang = 0.5 * q.ang;
    s = sin( hqang );
    
    retval.x = s * q.nx;
    retval.y = s * q.ny;
    retval.z = s * q.nz;
    retval.w = cos( hqang );
    
    /* q.qmag = hypotenuse4( q.x,  q.y,  q.z,  q.w  ); */
    return retval;
} // Quat convertRotToQuat( Quat q )

Quat raaToQuat( const Real raa[3], Real angle )
{
    Quat input;

    input.nx = raa[0];
    input.ny = raa[1];
    input.nz = raa[2];
    input.ang = angle;

    return convertRotToQuat( input );
} // Quat raaToQuat( Real raa[4] )

Quat uniformQuat( void )
    // Generate a uniformly-distributed random quaternion (UDQ)
{
    double x0, r1, r2, t1, t2;  // for uniformly distributed quaternion calculation
    Quat q;

    /*
    **  This should produce a uniformly distributed quaternion, according to
    **  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
    **  published by Academic Press, Inc., (1992)
    */
    t1 = genunf(0., TWOPI);
    // q.x = sin( t1 ) * (  r1 = ( (genunf(0., 1.) < 0.5) ?  (-1.) : (+1.) ) * sqrt( 1. - (x0 = genunf(0., 1.)) )  );  // random sign version
    q.x = sin( t1 ) * (  r1 = sqrt( 1. - (x0 = genunf(0., 1.)) )  );  // strict Shoemake version
    q.y = cos( t1 ) * r1;
    t2 = genunf(0., TWOPI);
    // q.z = sin( t2 ) * (  r2 = ( (genunf(0., 1.) < 0.5) ?  (-1.) : (+1.) ) * sqrt( x0 )  );  // random sign version
    q.z = sin( t2 ) * (  r2 = sqrt( x0 )  );  // strict Shoemake version
    q.w = cos( t2 ) * r2;

    return q;
}

Quat uniformQuatByAmount( Real amount )
    // returns a quaternion from a random axis and specified angle
    // amount is an angle in radians
{
    Quat q = uniformQuat();
    q = convertQuatToRot( q );
    q = axisRadianToQuat( q.nx, q.ny, q.nz, amount );
    return q;
}

void unitQuat2rotation( Quat *q )
    // Convert from a unit quaternion to a rotation about an unit 3D-vector
{
    double inv_sin_half_ang;

    q->ang = 2. * acos( q->w );
    inv_sin_half_ang = 1. / sin( 0.5 * q->ang );
    q->nx  = q->x * inv_sin_half_ang; 
    q->ny  = q->y * inv_sin_half_ang; 
    q->nz  = q->z * inv_sin_half_ang; 
    
    return;
}

void print_q_reorient_message( FILE *logFile, Quat q_reorient )
    // Print message about q_reorient
{
    pr( logFile, "\nRe-orienting the ligand using the following axis (nx, ny, nz) and angle values:\n");
    pr( logFile, "NEWDPF   reorient %.3lf %.3lf %.3lf %.2lf\n",
        q_reorient.nx, q_reorient.ny, q_reorient.nz, RadiansToDegrees( q_reorient.ang ) );

    pr( logFile, "\n");
    pr( logFile, "q_reorient:\n");
    printQuat( logFile, q_reorient );
    pr( logFile, "\n");

    return;
} // Print message about q_reorient

void create_random_orientation( Quat *ptr_quat ) 
{
    // Generate a random initial orientation for the ligand
    Quat q_random;
    // Generate a uniformly-distributed quaternion:
    // setting the x,y,z,w components
    q_random = uniformQuat();
    ptr_quat->x = q_random.x;
    ptr_quat->y = q_random.y;
    ptr_quat->z = q_random.z;
    ptr_quat->w = q_random.w;
    // Update the (nx,ny,nz,ang) components of the quaternion, ptr_quat:
    *ptr_quat = convertQuatToRot( *ptr_quat );
}

Quat conjugate( const Quat q )
{
    Quat conj;

    conj.x = -q.x;
    conj.y = -q.y;
    conj.z = -q.z;
    conj.w = q.w;

    return conj;
}

Quat inverse( const Quat q )
{
    register Quat conj, inv;
    register double inv_squared_magnitude;

    conj = conjugate( q );

    inv_squared_magnitude = 1. / sqhypotenuse4( conj.x, conj.y, conj.z, conj.w );

    inv.x = conj.x * inv_squared_magnitude;
    inv.y = conj.y * inv_squared_magnitude;
    inv.z = conj.z * inv_squared_magnitude;
    inv.w = conj.w * inv_squared_magnitude;

    return inv;
}

Quat slerp0( const Quat q1, const Quat q2, const double u )
    // See: Shoemake, K. (1985), "Animating Rotation with Quaternion Curves", 
    //      Computer Graphics, 19 (3): 245-254
    //
    // A formula for spherical linear interpolation from q1 to
    // q2, with parameter u moving from 0 to 1.
{
    Quat slerp;

    assert( u >= 0.  &&  u <= 1. );

    // q1 . q2 = cos( theta)
    double theta = acos( q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w );

    double w1 = sin( (1. - u) * theta ) / sin( theta );
    double w2 = sin(     u    * theta ) / sin( theta );

    slerp.x = w1 * q1.x  +  w2 * q2.x;
    slerp.y = w1 * q1.y  +  w2 * q2.y;
    slerp.z = w1 * q1.z  +  w2 * q2.z;
    slerp.w = w1 * q1.w  +  w2 * q2.w;

    return slerp;
}

Quat slerp1( const Quat qa, const Quat qb, const double t )
    // See Martin Baker's web site
    // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
{
    Quat qm; // quaternion to return
#ifdef ASSERTQUATOK
    assertQuatOK(qa);
    assertQuatOK(qb);
#endif

    assert( t >= 0.  &&  t <= 1. );

	// Calculate angle between them.
	double cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
	// if qa=qb or qa=-qb then theta = 0 and we can return qa
	if (fabs(cosHalfTheta) >= 1.0){
		qm.w = qa.w;qm.x = qa.x;qm.y = qa.y;qm.z = qa.z;
#ifdef ASSERTQUATOK
        assertQuatOK(qm);
#endif
		return qm;
	}
	// Calculate temporary values.
	double halfTheta = acos(cosHalfTheta);
	double sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
	// if theta = 180 degrees then result is not fully defined
	// we could rotate around any axis normal to qa or qb
	if (fabs(sinHalfTheta) < 0.001){ // fabs is floating point absolute
		qm.w = (qa.w * 0.5 + qb.w * 0.5);
		qm.x = (qa.x * 0.5 + qb.x * 0.5);
		qm.y = (qa.y * 0.5 + qb.y * 0.5);
		qm.z = (qa.z * 0.5 + qb.z * 0.5);
#ifdef DEBUG_MUTATION
        printQuat_q( logFile, qm );
        fprintf( logFile, "slerp:  WARNING!  theta = 180 degrees   " );
        printQuat_q( logFile, qm );
//        fflush(logFile);
#endif
#ifdef ASSERTQUATOK
        assertQuatOK(qm);
#endif
		return qm;
	}
	double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
	double ratioB = sin(t * halfTheta) / sinHalfTheta; 
	//calculate Quaternion.
	qm.w = qa.w * ratioA + qb.w * ratioB;
	qm.x = qa.x * ratioA + qb.x * ratioB;
	qm.y = qa.y * ratioA + qb.y * ratioB;
	qm.z = qa.z * ratioA + qb.z * ratioB;
#ifdef ASSERTQUATOK
    assertQuatOK(qm);
#endif
	return qm;
}

Quat slerp( const Quat qa, const Quat qb, const double t )
    // Adapted from code by John W. Ratcliff mailto:jratcliff@infiniplex.net
    // See http://codesuppository.blogspot.com/2006/03/matrix-vector-and-quaternion-library.html
/*  
** 
** Copyright (c) 2007 by John W. Ratcliff mailto:jratcliff@infiniplex.net
**
** The MIT license:
**
** Permission is hereby granted, free of charge, to any person obtaining a copy 
** of this software and associated documentation files (the "Software"), to deal 
** in the Software without restriction, including without limitation the rights 
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
** copies of the Software, and to permit persons to whom the Software is furnished 
** to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all 
** copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
** CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**
*/
{
    Quat qm; // quaternion to return
    Quat qb_local;
    double halfTheta, sinHalfTheta;
    double ratioA, ratioB;
#ifdef ASSERTQUATOK
    assertQuatOK(qa);
    assertQuatOK(qb);
#endif

    assert( t >= 0.  &&  t <= 1. );

	// Calculate angle between them.
	double cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
    // Ensure we choose the shorter angular displacement between qa and qb:
	if (cosHalfTheta < 0.) {
        cosHalfTheta = -cosHalfTheta;
        qb_local.w = -qb.w;
        qb_local.x = -qb.x;
        qb_local.y = -qb.y;
        qb_local.z = -qb.z;
	} else {
        qb_local = qb;
    }
#ifdef ASSERTQUATOK
    assertQuatOK(qb_local);
#endif
	// Calculate coefficients
    if ((1. - cosHalfTheta) > 1e-6) {
        // standard case (slerp)
        halfTheta = acos(cosHalfTheta);
        sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
        ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
        ratioB = sin(t * halfTheta) / sinHalfTheta; 
    } else {
        // qa and qb (qb_local) are very close, so we can do a linear interpolation
        ratioA = 1 - t ;
        ratioB = t; 
    }
	// Calculate final values
	qm.w = qa.w * ratioA + qb_local.w * ratioB;
	qm.x = qa.x * ratioA + qb_local.x * ratioB;
	qm.y = qa.y * ratioA + qb_local.y * ratioB;
	qm.z = qa.z * ratioA + qb_local.z * ratioB;
#ifdef ASSERTQUATOK
    assertQuatOK(qm);
#endif
	return qm;
}

Quat axisRadianToQuat( const Real ax, const Real ay, const Real az, const Real angle )
{
    Real raa[3] = { ax, ay, az };
    return raaToQuat( raa, angle );
}

Quat axisDegreeToQuat( const Real ax, const Real ay, const Real az, const Real angle )
{
    Real raa[3] = { ax, ay, az };
    return raaToQuat( raa, DegreesToRadians( angle ) );
}

Quat quatComponentsToQuat( const Real qx, const Real qy, const Real qz, const Real qw )
{
    Quat Q;
    Q.x = qx;
    Q.y = qy;
    Q.z = qz;
    Q.w = qw;
    return normQuat( Q );
}

const Quat identityQuat()
{
    Quat Q;
    Q.x = Q.y = Q.z = 0.;
    Q.w = 1.;
    return Q;
}

/* Radians */
#define ONE_ROTATION TWOPI // Degrees // #define ONE_ROTATION 360.
#define HALF_ROTATION PI // Degrees // #define HALF_ROTATION 180.

/* Angles that go from -half-a-rotation to half-a-rotation */
#define MIN_ANGLE -HALF_ROTATION // Angles that go from 0 to one-rotation // #define MIN_ANGLE 0.
#define MAX_ANGLE HALF_ROTATION // Angles that go from 0 to one-rotation // #define MAX_ANGLE ONE_ROTATION

Real a_range_reduction( Real a )
{
    if (a <= MIN_ANGLE) {
        do a += ONE_ROTATION;
        while (a <= MIN_ANGLE);
    } else if (a >= MAX_ANGLE) {
        do a -= ONE_ROTATION;
        while (a >= MAX_ANGLE);
    }
    return a;
}

Real alerp( Real a, Real b, Real fract )
{
    // if fract==0, return a
    // if fract==1, return b
    Real delta;
    a = a_range_reduction( a );
    b = a_range_reduction( b );
    delta = b - a;
    if (delta > HALF_ROTATION) {
        delta -= ONE_ROTATION;
    } else if (delta < -HALF_ROTATION) {
        delta += ONE_ROTATION;
    }
    return a_range_reduction( a + delta*fract );
}

/* test for alerp and a_range_reduction
int main() {
    Real start = -ONE_ROTATION;
    Real stop = ONE_ROTATION;
    Real step = ONE_ROTATION/8.;
    Real i, j;
    for (i=start; i<stop; i=i+step) {
        printf(" %.3f:\n", i); 
        for (j=start; j<stop; j=j+step) {
            printf("      %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", j, alerp(i,j,0.0), alerp(i,j,0.1), alerp(i,j,0.5), alerp(i,j,0.9), alerp(i,j,1.0));
        }
    }
    return 0;
}
*/

/* EOF */
