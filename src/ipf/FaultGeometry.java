/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Methods for fault geometry in seismic image coordinates.
 * <p>
 * Vectors u, v, and w denote fault dip, strike, and normal vectors,
 * respectively. These vectors form a right-handed coordinate system
 * in which u = v x w, v = w x u, and w = u x v, where "x" denotes 
 * the vector cross product.
 * <p>
 * Likewise, strike and dip angles are defined by the right-hand-rule. That
 * is, if fingers of the right hand are rotated to lie in the fault plane and
 * to point downward in the fault dip direction, then the thumb of that hand
 * points in the fault strike direction. Strike and dip angles are measured in
 * degrees.
 * <p>
 * The components of all vectors, (e.g., u1, u2, and u3) correspond to the
 * 1st, 2nd, and 3rd dimensions of an image. The 1st image axis is vertical,
 * and the 2nd and 3rd axes typically correspond to inline and crossline
 * seismic survey directions. If e1, e2, and e3 are unit vectors aligned
 * with the image axes, then e1 = e3 x e2, e2 = e1 x e3, and e3 = e2 x e1.
 * <p>
 * The fault strike angle phi is measured clockwise in a horizontal plane from
 * image axis x3 to image axis x2. Likewise, the fault dip angle theta is
 * positive and increases downward from the horizontal plane. Note that the
 * true (conventional) strike angle can be easily computed from the image
 * strike angle, if the orientation of the seismic survey is known.
 * <p>
 * Methods with strike angle parameters will accept any value. Methods that
 * return strike angles will return values in the range [0,360] degrees.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.29
 */
public class FaultGeometry {

  /**
   * Returns fault dip vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {u1,u2,u3} of components for dip vector.
   */
  public static float[] faultDipVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float u1 = (float)( st);
    float u2 = (float)( ct*cp);
    float u3 = (float)(-ct*sp);
    return new float[]{u1,u2,u3};
  }

  /**
   * Returns fault strike vector for specified strike and dip angles.
   * The dip angle theta is not used, but is provided for consistency.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {v1,v2,v3} of components for strike vector.
   */
  public static float[] faultStrikeVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double cp = cos(p);
    double sp = sin(p);
    float v1 = 0.0f;
    float v2 = (float)sp;
    float v3 = (float)cp;
    return new float[]{v1,v2,v3};
  }

  /**
   * Returns fault normal vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {w1,w2,w3} of components for normal vector.
   */
  public static float[] faultNormalVectorFromStrikeAndDip(
      double phi, double theta) {
    double p = toRadians(phi);
    double t = toRadians(theta);
    double cp = cos(p);
    double sp = sin(p);
    double ct = cos(t);
    double st = sin(t);
    float w1 = (float)(-ct);
    float w2 = (float)( st*cp);
    float w3 = (float)(-st*sp);
    return new float[]{w1,w2,w3};
  }

  /**
   * Returns fault strike angle for specified fault dip vector.
   * The components u2 and u3 must not both be zero; that is, the
   * fault dip vector cannot be vertical.
   * @param u1 1st component of fault dip vector.
   * @param u2 2nd component of fault dip vector.
   * @param u3 3rd component of fault dip vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromDipVector(float u1, float u2, float u3) {
    Check.argument(u2!=0.0f || u3!=0.0f,"dip vector is not vertical");
    return range360(toDegrees(atan2(-u3,u2)));
  }

  /**
   * Returns fault strike angle for specified fault dip vector.
   * The components u2 and u3 must not both be zero; that is, the
   * fault dip vector cannot be vertical.
   * @param u array {u1,u2,u3} of components of fault dip vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromDipVector(float[] u) {
    return faultStrikeFromDipVector(u[0],u[1],u[2]);
  }

  /**
   * Returns fault dip angle for specified fault dip vector.
   * @param u1 1st component of fault dip vector.
   * @param u2 2nd component of fault dip vector.
   * @param u3 3rd component of fault dip vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromDipVector(float u1, float u2, float u3) {
    return toDegrees(asin(u1));
  }

  /**
   * Returns fault dip angle for specified fault dip vector.
   * @param u array {u1,u2,u3} of components of fault dip vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromDipVector(float[] u) {
    return faultDipFromDipVector(u[0],u[1],u[2]);
  }

  /**
   * Returns fault strike angle for specified fault strike vector.
   * @param v1 1st component of fault strike vector.
   * @param v2 2nd component of fault strike vector.
   * @param v3 3rd component of fault strike vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromStrikeVector(
      float v1, float v2, float v3) {
    return range360(toDegrees(atan2(v2,v3)));
  }

  /**
   * Returns fault strike angle for specified fault strike vector.
   * @param v array {v1,v2,v3} of components of fault strike vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromStrikeVector(float[] v) {
    return faultStrikeFromStrikeVector(v[0],v[1],v[2]);
  }

  /**
   * Returns fault strike angle for specified fault normal vector.
   * The components w2 and w3 must not both be zero; that is, the
   * fault plane cannot be horizontal.
   * @param w1 1st component of fault normal vector.
   * @param w2 2nd component of fault normal vector.
   * @param w3 3rd component of fault normal vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromNormalVector(
      float w1, float w2, float w3) {
    Check.argument(w2!=0.0f || w3!=0.0f,"normal vector is not vertical");
    return range360(toDegrees(atan2(-w3,w2)));
  }

  /**
   * Returns fault strike angle for specified fault normal vector.
   * The components w2 and w3 must not both be zero; that is, the
   * fault plane cannot be horizontal.
   * @param w array {w1,w2,w3} of components of fault normal vector.
   * @return fault strike angle, in degrees.
   */
  public static float faultStrikeFromNormalVector(float[] w) {
    return faultStrikeFromNormalVector(w[0],w[1],w[2]);
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param w1 1st component of fault normal vector.
   * @param w2 2nd component of fault normal vector.
   * @param w3 3rd component of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormalVector(float w1, float w2, float w3) {
    return toDegrees(acos(-w1));
  }

  /**
   * Returns fault dip angle for specified fault normal vector.
   * @param w array {w1,w2,w3} of components of fault normal vector.
   * @return fault dip angle, in degrees.
   */
  public static float faultDipFromNormalVector(float[] w) {
    return faultDipFromNormalVector(w[0],w[1],w[2]);
  }

  /**
   * Returns the cross-product of two specified vectors.
   * @param u the 1st vector.
   * @param v the 2nd vector.
   * @return array {w1,w2,w3} of components for vector w = u x v.
   */
  public static float[] crossProduct(float[] u, float[] v) {
    float u1 = u[0], u2 = u[1], u3 = u[2];
    float v1 = v[0], v2 = v[1], v3 = v[2];
    float w1 = u3*v2-u2*v3;
    float w2 = u1*v3-u3*v1;
    float w3 = u2*v1-u1*v2;
    return new float[]{w1,w2,w3};
  }

  /**
   * Returns angle in range [0,360] degrees.
   * @param phi angle, in degrees.
   * @return angle in range [0,360] degrees.
   */
  public static float range360(double phi) {
    while (phi<0.0)
      phi += 360.0;
    while (phi>=360.0)
      phi -= 360.0;
    return (float)phi;
  }

  /**
   * Returns angle in range [-180,180] degrees.
   * @param phi angle.
   * @return angle in range [-180,180] degrees.
   */
  public static float range180(double phi) {
    while (phi<-180.0)
      phi += 360.0;
    while (phi>180.0)
      phi -= 360.0;
    return (float)phi;
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    float[] fp = {  0.0f, 90.0f,180.0f,270.0f,
                    0.0f, 90.0f,180.0f,270.0f};
    float[] ft = { 90.0f, 90.0f, 90.0f, 90.0f,
                   89.0f, 89.0f, 89.0f, 89.0f};
    for (int i=0; i<fp.length; ++i)
      test(fp[i],ft[i]);
  }
  public static void test(float phia, float thetaa) {
    float[] ua = faultNormalVectorFromStrikeAndDip(phia,thetaa);
    float phib = faultStrikeFromNormalVector(ua);
    float thetab = faultDipFromNormalVector(ua);
    float[] ub = faultNormalVectorFromStrikeAndDip(phib,thetab);
    assertEqual(ua,ub);
  }
  public static void assertEqual(float x, float y) {
    assert abs(x-y)<0.01f;
  }
  public static void assertEqual(float[] x, float[] y) {
    for (int i=0; i<x.length; ++i)
      assertEqual(x[i],y[i]);
  }
  public static void trace(String s) {
    System.out.println(s);
  }
}
