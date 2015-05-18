/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

/**
 * An abstract vector.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public interface Vec {

  /**
   * Returns the machine epsilon for this vector. Machine epsilon is
   * a measure of the precision with which vector elements are stored. 
   * It is approximately equal to the smallest positive number such 
   * that 1!=1+epsilon.
   * @return the machine epsilon.
   */
  public double epsilon();

  /**
   * Returns a clone of this vector.
   */
  public Vec clone();

  /**
   * Returns the dot product of this vector with that vector.
   * @param vthat that vector.
   * @return the dot product.
   */
  public double dot(Vec vthat);

  /**
   * Returns the L2 norm of this vector.
   * The L2 norm is also called the Euclidean norm. It equals the
   * square root of the sum of squared elements of this vector.
   * @return the L2 norm.
   */
  public double norm2();

  /**
   * Zeros all elements of this vector.
   */
  public void zero();

  /**
   * Scales this vector by the specified factor.
   * @param s the scale factor.
   */
  public void scale(double s);

  /**
   * Updates this vector by computing vthis = vthis*sthis + vthat*sthat.
   * @param sthis factor by which to scale this vector.
   * @param vthat that vector.
   * @param sthat factor by which to scale that vector.
   */
  public void add(double sthis, Vec vthat, double sthat);
}
