/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Parallel.*;

/**
 * A vector represented by a 3D array[n3][n2][n1] of floats.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.29
 */
public class VecArrayFloat3 implements Vec {

  /**
   * Constructs a zero vector with specified dimensions.
   * @param n1 the number of floats in the 1st dimension.
   * @param n2 the number of floats in the 2nd dimension.
   * @param n3 the number of floats in the 3rd dimension.
   */
  public VecArrayFloat3(int n1, int n2, int n3) {
    _a = new float[n3][n2][n1];
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
  }

  /**
   * Constructs a vector that wraps the specified array of floats.
   * @param a the array of floats; by reference, not by copy.
   */
  public VecArrayFloat3(float[][][] a) {
    _a = a;
    _n1 = a[0][0].length;
    _n2 = a[0].length;
    _n3 = a.length;
  }

  /**
   * Gets the array of floats wrapped by this vector.
   * @return the array of floats; by reference, not by copy.
   */
  public float[][][] getArray() {
    return _a;
  }

  /**
   * Gets the number of floats in the 1st array dimension.
   * @return the number of floats in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of floats in the 2nd array dimension.
   * @return the number of floats in the 2nd dimension.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the number of floats in the 3rd array dimension.
   * @return the number of floats in the 3rd dimension.
   */
  public int getN3() {
    return _n3;
  }

  public double epsilon() {
    return Math.ulp(1.0f);
  }

  public VecArrayFloat3 clone() {
    VecArrayFloat3 v = new VecArrayFloat3(_n1,_n2,_n3);
    scopy(_a,v._a);
    return v;
  }

  public double dot(Vec vthat) {
    float[][][] athis = _a;
    float[][][] athat = ((VecArrayFloat3)vthat)._a;
    return sdot(athis,athat);
  }

  public double norm2() {
    return Math.sqrt(sdot(_a,_a));
  }

  public void zero() {
    szero(_a);
  }

  public void scale(double s) {
    sscal((float)s,_a);
  }

  public void add(double sthis, Vec vthat, double sthat) {
    float fthis = (float)sthis;
    float fthat = (float)sthat;
    float[][][] athis = _a;
    float[][][] athat = ((VecArrayFloat3)vthat)._a;
    if (fthis==1.0f) {
      saxpy(fthat,athat,athis);
    } else if (fthat==1.0f) {
      sxpay(fthis,athat,athis);
    } else {
      saxpby(fthat,athat,fthis,athis);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _a;
  private int _n1,_n2,_n3;

  // Zeros array x.
  private static void szero(float[] x) {
    ArrayMath.zero(x);
  }
  private static void szero(float[][] x) {
    ArrayMath.zero(x);
  }
  private void szero(final float[][][] x) {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      szero(x[i3]);
    }});
  }

  // Copys array x to array y.
  private void scopy(float[] x, float[] y) {
    ArrayMath.copy(x,y);
  }
  private void scopy(float[][] x, float[][] y) {
    ArrayMath.copy(x,y);
  }
  private void scopy(final float[][][] x, final float[][][] y) {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      scopy(x[i3],y[i3]);
    }});
  }

  // Returns the dot product x'y.
  private double sdot(float[] x, float[] y) {
    int n1 = x.length;
    double d = 0.0;
    for (int i1=0; i1<n1; ++i1)
      d += x[i1]*y[i1];
    return d;
  }
  private double sdot(float[][] x, float[][] y) {
    int n2 = x.length;
    double d = 0.0;
    for (int i2=0; i2<n2; ++i2)
      d += sdot(x[i2],y[i2]);
    return d;
  }
  private double sdot(final float[][][] x, final float[][][] y) {
    int n3 = x.length;
    double d = reduce(n3,new ReduceInt<Double>() {
      public Double compute(int i3) {
        return sdot(x[i3],y[i3]);
      }
      public Double combine(Double a, Double b) {
        return a+b;
      }
    });
    return d;
  }

  // Computes x = a*x.
  private void sscal(float a, float[] x) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      x[i1] *= a;
  }
  private void sscal(float a, float[][] x) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      sscal(a,x[i2]);
  }
  private void sscal(final float a, final float[][][] x) {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      sscal(a,x[i3]);
    }});
  }

  // Computes y = y + a*x.
  private void saxpy(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += a*x[i1];
  }
  private void saxpy(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      saxpy(a,x[i2],y[i2]);
  }
  private void saxpy(
    final float a, final float[][][] x, final float[][][] y)
  {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      saxpy(a,x[i3],y[i3]);
    }});
  }

  // Computes y = x + a*y.
  private void sxpay(float a, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*y[i1]+x[i1];
  }
  private void sxpay(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      sxpay(a,x[i2],y[i2]);
  }
  private void sxpay(
    final float a, final float[][][] x, final float[][][] y)
  {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      sxpay(a,x[i3],y[i3]);
    }});
  }

  // Computes y = a*x + b*y.
  private void saxpby(float a, float[] x, float b, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] = a*x[i1]+b*y[i1];
  }
  private void saxpby(float a, float[][] x, float b, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2)
      saxpby(a,x[i2],b,y[i2]);
  }
  private void saxpby(
    final float a, final float[][][] x, final float b, final float[][][] y)
  {
    int n3 = x.length;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      saxpby(a,x[i3],b,y[i3]);
    }});
  }
}
