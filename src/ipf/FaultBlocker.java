/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import java.util.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static ipf.FaultGeometry.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.09.16
 */
public class FaultBlocker {

  /**
   * Returns fault blocks for specified fault images.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes and dips.
   * @return array of fault blocks.
   */
  public float[][][] findBlocks(float[][][][] flpt) {
    return blocks(flpt);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Returns fault blocks.
  private static float[][][] blocks(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Compute right-hand-side.
    float[][][] r = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int m3 = (i3>0)?i3-1:0;
      for (int i2=0; i2<n2; ++i2) {
        int m2 = (i2>0)?i2-1:0;
        for (int i1=0; i1<n1; ++i1) {
          int m1 = (i1>0)?i1-1:0;
          float fl = f[i3][i2][i1];
          float fp = p[i3][i2][i1];
          float ft = t[i3][i2][i1];
          float[] fn = faultNormalVectorFromStrikeAndDip(fp,ft);
          float fn1 = fn[0];
          float fn2 = fn[1];
          float fn3 = fn[2];
          float fl1 = 0.5f*(f[i3][i2][i1]+f[i3][i2][m1]);
          float fl2 = 0.5f*(f[i3][i2][i1]+f[i3][m2][i1]);
          float fl3 = 0.5f*(f[i3][i2][i1]+f[m3][i2][i1]);
          float f1 = fl1*fn1;
          float f2 = fl2*fn2;
          float f3 = fl3*fn3;
          r[i3][i2][i1] += f1+f2+f3;
          r[i3][i2][m1] -= f1;
          r[i3][m2][i1] -= f2;
          r[m3][i2][i1] -= f3;
        }
      }
    }

    // Solve for fault blocks.
    A3 a = new A3();
    CgSolver cs = new CgSolver(0.01,100);
    float[][][] b = new float[n3][n2][n1];
    VecArrayFloat3 vb = new VecArrayFloat3(b);
    VecArrayFloat3 vr = new VecArrayFloat3(r);
    cs.solve(a,vr,vb);
    return b;
  }
  private static class A3 implements CgSolver.A {
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat3 v3x = (VecArrayFloat3)vx;
      VecArrayFloat3 v3y = (VecArrayFloat3)vy;
      float[][][] x = v3x.getArray();
      float[][][] y = v3y.getArray();
      zero(y);
      _ldk.apply(x,y);
    }
    private LocalDiffusionKernel _ldk =
      new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21);
  }
}
