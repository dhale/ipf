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
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;
import static ipf.FaultGeometry.*;

/**
 * Uses image samples alongside fault skins to estimate fault dip-slips.
 * <p>
 * Before estimating slips, the seismic images should be smoothed along
 * reflectors, but not across faults. (Thinned fault likelihoods can be used
 * to stop smoothing at potential fault locations.) This smoothing does what
 * seismic interpreters do visually when estimating fault throws. In effect,
 * such smoothing enables us to use samples of a seismic image located away
 * from faults to estimate slips along those faults.
 * <p>
 * However, even after this smoothing, we may not want to use image samples
 * located immediately adjacent to faults, as these samples may exhibit poor
 * estimates of reflector slopes and smoothing filter artifacts. Therefore,
 * image samples at some small horizontal distance (e.g., two samples away)
 * from fault skins are used to estimate slips. Specified reflector slopes are
 * then used to account for this small offset when computing fault slips.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.05
 */
public class FaultSlipper {

  /**
   * Constructs a fault slipper for the specified seismic image and slopes.
   * @param gs seismic image, smoothed up to (but not across) faults.
   * @param p2 slopes in the 2nd dimension.
   * @param p3 slopes in the 3rd dimension.
   */
  public FaultSlipper(float[][][] gs, float[][][] p2, float[][][] p3) {
    _gs = gs;
    _p2 = p2;
    _p3 = p3;
    _n1 = gs[0][0].length;
    _n2 = gs[0].length;
    _n3 = gs.length;
  }

  /**
   * Sets the offset distance used when estimating slips.
   * The default offset is 2 samples.
   * @param offset the offset, in samples.
   */
  public void setOffset(double offset) {
    _offset = (float)abs(offset);
  }

  /**
   * Enables or disables a zero-slope reflector assumption. For educational
   * use, only. Ignoring reflector slopes can yield significant errors in
   * estimated dip slips. The default is false.
   * @param zeroSlope true, to assume zero slopes; false, otherwise.
   */
  public void setZeroSlope(boolean zeroSlope) {
    _zeroSlope = zeroSlope;
  }

  /**
   * Computes fault dip slips for the specified skins.
   * Specified bounds on throw are used to reduce unnecessary computation.
   * These bounds need not be precise, and are not enforced. If the fault
   * skins represent reverse faults (with negative throws), the specified
   * estimate for minimum fault throw should be negative.
   * @param skins array of skins for which to compute dip slips.
   * @param smin an estimate for minimum fault throw, in samples.
   * @param smax an estimate for maximum fault throw, in samples.
   */
  public void computeDipSlips(FaultSkin[] skins, double smin, double smax) {
    for (FaultSkin skin:skins)
      computeDipSlips(skin,smin,smax);
  }

  /**
   * Computes fault dip slips for all cells in the specified skin.
   * Specified bounds on throw are used to reduce unnecessary computation.
   * These bounds need not be precise, and are not enforced. If the fault
   * skins represent reverse faults (with negative throws), the specified
   * estimate for minimum fault throw should be negative.
   * @param skin the skin for which to compute shifts.
   * @param smin an estimate for minimum fault throw, in samples.
   * @param smax an estimate for maximum fault throw, in samples.
   */
  public void computeDipSlips(FaultSkin skin, double smin, double smax) {
    Check.argument(smax>=0.0f,"smax not less than zero");
    FaultCell[][] cab = skin.getCellsAB();
    FaultCell[][] clr = skin.getCellsLR();
    int lmin = round((float)smin-2*_offset); // because shifts != throws
    int lmax = round((float)smax+2*_offset); // TODO: use pmax*_offset?
    DynamicWarping dw = new DynamicWarping(lmin,lmax);
    dw.setStrainMax(0.25,0.25); // TODO: always 0.25? goes with 4 below?
    computeAlignmentErrors(skin,lmin,lmax,_offset,_gs);
    extrapolateAlignmentErrors(lmin,lmax,cab);
    computeShifts(dw,cab,clr);
    clearErrors(skin);
    for (int nsmooth=0; nsmooth<2; ++nsmooth) // TODO: 2?
      smoothShifts(skin);
    computeDipSlips(skin);
  }

  /**
   * Gets images of fault dip-slip vectors for an array of fault skins.
   * Assumes that dip-slip vectors have already been computed and stored
   * in all skin cells.
   * <p>
   * Dip-slip vectors are non-zero for only the hanging wall (plus side)
   * of faults. Slips on the footwall (minus side) are zero. Image samples
   * not adjacent to faults are marked using a specified value.
   * @param skins array of fault skins.
   * @param smark the mark for slips not adjacent to a fault.
   * @return array {s1,s2,s3} of components of dip slips.
   */
  public float[][][][] getDipSlips(FaultSkin[] skins, float smark) {
    int n1 = _n1, n2 = _n2, n3 = _n3;
    float[][][] s1 = new float[n3][n2][n1];
    float[][][] s2 = new float[n3][n2][n1];
    float[][][] s3 = new float[n3][n2][n1];
    float[][][] ss = new float[n3][n2][n1];

    // Initially set all slip vectors to the specified mark.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          s1[i3][i2][i1] = smark;
          s2[i3][i2][i1] = smark;
          s3[i3][i2][i1] = smark;
        }
      }
    }

    // For all cells in all skins, ...
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {

        // Get sample indices for the minus and plus sides of the cell.
        int i1 = cell.i1;
        int i2m = cell.i2m;
        int i3m = cell.i3m;
        int i2p = cell.i2p;
        int i3p = cell.i3p;

        // If slip on the minus side has not been set, zero it.
        if (s1[i3m][i2m][i1]==smark) {
          s1[i3m][i2m][i1] = 0.0f;
          s2[i3m][i2m][i1] = 0.0f;
          s3[i3m][i2m][i1] = 0.0f;
        }

        // Set or accumulate slip on the plus side.
        if (s1[i3p][i2p][i1]==smark) {
          s1[i3p][i2p][i1]  = cell.s1;
          s2[i3p][i2p][i1]  = cell.s2;
          s3[i3p][i2p][i1]  = cell.s3;
          ss[i3p][i2p][i1]  = 1.0f;
        } else {
          s1[i3p][i2p][i1] += cell.s1;
          s2[i3p][i2p][i1] += cell.s2;
          s3[i3p][i2p][i1] += cell.s3;
          ss[i3p][i2p][i1] += 1.0f;
        }
      }
    }

    // Where more than one slip was accumulated, compute the average.
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if (ss[i3][i2][i1]>1.0f) {
            float si = 1.0f/ss[i3][i2][i1];
            s1[i3][i2][i1] *= si;
            s2[i3][i2][i1] *= si;
            s3[i3][i2][i1] *= si;
          }
        }
      }
    }
    return new float[][][][]{s1,s2,s3};
  }

  /**
   * Interpolates specified dip-slip vectors.
   * @param s array {s1,s2,s3} of dip-slip vectors.
   * @param smark the mark for slips not adjacent to a fault.
   * @return interpolated dip-slip vectors.
   */
  public static float[][][][] interpolateDipSlips(
      float[][][][] s, float smark) {
    int n1 = s[0][0][0].length;
    int n2 = s[0][0].length;
    int n3 = s[0].length;
    float[][][] sp = new float[n3][n2][n1];
    float[][][] st = new float[n3][n2][n1];
    short[][][] k1 = new short[n3][n2][n1];
    short[][][] k2 = new short[n3][n2][n1];
    short[][][] k3 = new short[n3][n2][n1];
    float[][][][] sq = new float[3][n3][n2][n1];
    ClosestPointTransform cpt = new ClosestPointTransform();
    cpt.apply(smark,s[0],st,k1,k2,k3);
    clip(0.0f,100.0f,st,st);
    LocalDiffusionKernel.Stencil stencil = LocalDiffusionKernel.Stencil.D21;
    LocalDiffusionKernel ldk = new LocalDiffusionKernel(stencil);
    BlendedGridder3 bg = new BlendedGridder3();
    bg.setBlendingKernel(ldk);
    bg.setSmoothness(0.5);
    for (int is=0; is<3; ++is) {
      float[][][] si = s[is];
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            int j1 = k1[i3][i2][i1];
            int j2 = k2[i3][i2][i1];
            int j3 = k3[i3][i2][i1];
            sp[i3][i2][i1] = si[j3][j2][j1];
          }
        }
      }
      bg.gridBlended(st,sp,sq[is]);
    }
    return sq;
  }

  /**
   * Unfaults an image using interpolated dip-slip vectors.
   * @param s array {s1,s2,s3} of interpolated dip-slip vectors.
   * @param g image to be unfaulted.
   * @return unfaulted image.
   */
  public static float[][][] unfault(float[][][][] s, final float[][][] g) {
    final int n1 = g[0][0].length;
    final int n2 = g[0].length;
    final int n3 = g.length;
    final float[][][] s1 = s[0];
    final float[][][] s2 = s[1];
    final float[][][] s3 = s[2];
    final float[][][] gs = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1 = i1+s1[i3][i2][i1];
          float x2 = i2+s2[i3][i2][i1];
          float x3 = i3+s3[i3][i2][i1];
          gs[i3][i2][i1] = si.interpolate(
              n1,1.0,0.0,
              n2,1.0,0.0,
              n3,1.0,0.0,
              g,x1,x2,x3);
        }
      }
    }});
    return gs;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float[][][] _gs,_p2,_p3; // seismic image (smoothed) and slopes
  private int _n1,_n2,_n3; // image dimensions
  private float _offset = 2.0f; // horizontal offset (distance to faults)
  private boolean _zeroSlope; // if true, assume reflectors have zero slope

  /**
   * Computes alignment errors and initializes shifts for specified skin.
   */
  private static void computeAlignmentErrors(
      FaultSkin skin, int lmin, int lmax, float offset, float[][][] f) {
    for (FaultCell cell:skin)
      computeAlignmentErrors(cell,lmin,lmax,offset,f);
  }

  /**
   * Computes minus-plus alignment errors for one cell. These errors
   * correspond to differences between the sample value on the minus side of
   * the cell and those for the plus sides of cells up and down dip from the
   * cell.
   * <p> 
   * For lags where image sample values are unavailable (e.g., near image
   * boundaries), errors are extrapolated from other lags, but are negated, so
   * that extrapolated errors can be detected and modified later, after errors
   * for all relevant cells have been computed.
   */
  private static void computeAlignmentErrors(
      FaultCell cell, int lmin, int lmax, float offset, float[][][] f) {
    Check.argument(lmin<=0,"lmin<=0");
    Check.argument(lmax>=0,"lmax>=0");
    int n1 = f[0][0].length;
    float[] y = new float[3];

    // New arrays for alignment errors.
    int lag0 = -lmin;
    float[] emp = cell.emp = new float[1+lmax-lmin];

    // Errors for lag zero.
    float d2 =  offset*cell.v3;
    float d3 = -offset*cell.v2;
    float y1 = cell.x1, y2 = cell.x2, y3 = cell.x3;
    float fm = imageValueAt(y1,y2-d2,y3-d3,f);
    float gp = imageValueAt(y1,y2+d2,y3+d3,f);
    float empl = emp[lag0] = alignmentError(fm,gp);

    // Errors for samples above; make any extrapolated errors negative.
    FaultCell ca = cell;
    int nlaga = min(-lmin,ca.i1);
    y1 = cell.x1; y2 = cell.x2; y3 = cell.x3;
    for (int ilag=1; ilag<=-lmin; ++ilag) {
      if (ilag<=nlaga) {
        y[0] = y1; y[1] = y2; y[2] = y3;
        ca = ca.walkUpDipFrom(y);
        y1 = y[0]; y2 = y[1]; y3 = y[2];
        d2 =  offset*ca.v3;
        d3 = -offset*ca.v2;
        gp = imageValueAt(y1,y2+d2,y3+d3,f);
        empl = emp[lag0-ilag] = alignmentError(fm,gp);
      } else {
        emp[lag0-ilag] = -empl;
      }
    }

    // Errors for samples below; make any extrapolated errors negative.
    FaultCell cb = cell;
    int nlagb = min(lmax,n1-1-cb.i1);
    y1 = cell.x1; y2 = cell.x2; y3 = cell.x3;
    for (int ilag=1; ilag<=lmax; ++ilag) {
      if (ilag<=nlagb) {
        y[0] = y1; y[1] = y2; y[2] = y3;
        cb = cb.walkDownDipFrom(y);
        y1 = y[0]; y2 = y[1]; y3 = y[2];
        d2 =  offset*cb.v3;
        d3 = -offset*cb.v2;
        gp = imageValueAt(y1,y2+d2,y3+d3,f);
        empl = emp[lag0+ilag] = alignmentError(fm,gp);
      } else {
        emp[lag0+ilag] = -empl;
      }
    }
  }

  /**
   * Extrapolates alignment errors emp where not computed. Errors that could
   * not be computed are negative, and are copies of errors for smaller lags
   * that could be computed. (Errors for zero lag can always be computed.) 
   * <p>
   * For each lag with a negative error, this method first attempts to
   * extrapolate using other errors for the same lag stored in cell nabors
   * above or below. This first extrapolation, where feasible, works well when
   * shifts vary slowly with depth.
   * <p> 
   * If this first extrapolation is infeasible, because the number of above
   * and below nabors for some lag is too small, then errors are extrapolated
   * using the errors already computed for other lags. Those errors are
   * already stored in the cells, but are negative, so in this second
   * extrapolation we simply change their sign.
   */
  private static void extrapolateAlignmentErrors(
      int lmin, int lmax, FaultCell[][] cab) {
    int nab = cab.length;

    // For all arrays of cells linked above-below, ...
    for (int iab=0; iab<nab; ++iab) {
      int mab = cab[iab].length;
      float[][] emp = new float[mab][];

      // Make arrays of errors for all cells in one above-below array.
      for (int jab=0; jab<mab; ++jab) {
        FaultCell c = cab[iab][jab];
        emp[jab] = c.emp;
      }

      // For each cell (each array of errors), ...
      for (int jab=0; jab<mab; ++jab) {

        // For all lags, ...
        for (int lag=lmin,ilag=0; lag<=lmax; ++lag,++ilag) {

          // The error for one cell and one lag.
          float empi = emp[jab][ilag];

          // If negative, this error was extrapolating using errors for
          // other lags, so search for an error with the same lag.
          if (empi<0.0f) {

            // If lag is negative, search below; otherwise, search above.
            if (lag<0) {
              for (int kab=jab; kab<mab && empi<0.0f; ++kab)
                empi = emp[kab][ilag];
            } else if (lag>0) {
              for (int kab=jab; kab>=0 && empi<0.0f; --kab)
                empi = emp[kab][ilag];
            }

            // If no good error found, use what we have (but made positive).
            if (empi<0.0f) 
              empi = -emp[jab][ilag];
          }

          // Update the error stored in the cell for one lag.
          emp[jab][ilag] = empi;
        } 
      }
    }
  }

  /**
   * Uses dynamic warping to compute minus-plus shifts. Assumes that
   * minus-plus errors have already been computed and stored in the cells
   * referenced in the specified above-below and left-right arrays.
   */
  private static void computeShifts(
      DynamicWarping dw, FaultCell[][] cab, FaultCell[][] clr) {
    
    // Arrays of arrays of errors, linked above and below.
    int nab = cab.length;
    float[][][] eab = new float[nab][][];
    for (int iab=0; iab<nab; ++iab) {
      int mab = cab[iab].length;
      eab[iab] = new float[mab][];
      for (int jab=0; jab<mab; ++jab) {
        FaultCell c = cab[iab][jab];
        eab[iab][jab] = c.emp;
      }
    }

    // Arrays of arrays of errors, linked left and right.
    int nlr = clr.length;
    float[][][] elr = new float[nlr][][];
    for (int ilr=0; ilr<nlr; ++ilr) {
      int mlr = clr[ilr].length;
      elr[ilr] = new float[mlr][];
      for (int jlr=0; jlr<mlr; ++jlr) {
        FaultCell c = clr[ilr][jlr];
        elr[ilr][jlr] = c.emp;
      }
    }


    // Smooth alignment errors in above-below and left-right directions.
    for (int ismooth=0; ismooth<2; ++ismooth) { // TODO: how many?
      dw.smoothErrors1(eab,eab);
      normalizeErrors(eab); // TODO: helpful?
      dw.smoothErrors1(elr,elr);
      normalizeErrors(elr); // TODO: helpful?
    }

    // Find shifts by accumulating once more and then backtracking.
    for (int iab=0; iab<nab; ++iab) {
      float[][] dab = dw.accumulateForward(eab[iab]);
      float[] s = dw.backtrackReverse(dab,eab[iab]);
      int mab = s.length;
      for (int jab=0; jab<mab; ++jab) {
        FaultCell c = cab[iab][jab];
        c.smp = s[jab];
      }
    }
  }
  
  /**
   * Normalizes errors to account for varying lengths of arrays of cells.
   * This normalization is different from that performed by dynamic warping
   * when smoothing alignment errors, because that normalization assumes that
   * all arrays in the array of arrays of alignment errors have the same
   * length. This assumption is false for the arrays cellsAB and cellsLR in
   * fault skins.
   */
  private static void normalizeErrors(float[][][] e) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2) {
      int n1 = e[i2].length;
      float scale = 1.0f/n1;
      for (int i1=0; i1<n1; ++i1) {
        float[] ei = e[i2][i1];
        int nlag = ei.length;
        for (int ilag=0; ilag<nlag; ++ilag)
          ei[ilag] *= scale;
      }
    }
  }

  /**
   * Smooths shifts in the specified skin.
   */
  private static void smoothShifts(FaultSkin skin) {
    FaultCell.Get1 getter = new FaultCell.Get1() { 
      public float get(FaultCell cell) { return cell.smp; }
    };
    FaultCell.Set1 setter = new FaultCell.Set1() { 
      public void set(FaultCell cell, float smp) { cell.smp = smp; }
    };
    skin.smooth1(getter,setter);
  }

  /**
   * Computes dip-slip vectors from vertical shifts for specified skin.
   */
  private void computeDipSlips(FaultSkin skin) {

    // For all cells in the skin, ...
    for (FaultCell cell:skin) {
      FaultCell cellBegin = cell;

      // Cell coordinates, dip vector, and shift.
      float x1 = cell.x1;
      float x2 = cell.x2;
      float x3 = cell.x3;
      float u1 = cell.u1;
      float u2 = cell.u2;
      float u3 = cell.u3;
      float smp = cell.smp;

      // Offset vector d.
      float d1 = 0.0f;
      float d2 =  _offset*cell.v3;
      float d3 = -_offset*cell.v2;

      // Reflector slopes at point x-d.
      float p2 = imageValueAt(x1-d1,x2-d2,x3-d3,_p2);
      float p3 = imageValueAt(x1-d1,x2-d2,x3-d3,_p3);

      // Unit-vector a normal to reflector.
      float a1 = 1.0f/sqrt(1.0f+p2*p2+p3*p3);
      float a2 = -p2*a1;
      float a3 = -p3*a1;

      // Slip adjustment for reflector slope on minus side.
      float am = (a1*d1+a2*d2+a3*d3)/(a1*u1+a2*u2+a3*u3);
      float u1m = am*u1;
      float u2m = am*u2;
      float u3m = am*u3;

      // Begin at cell location.
      float[] y = {cell.x1,cell.x2,cell.x3};

      // Walk down-dip (for normal fault) or up-dip (for reverse fault).
      if (smp>0.0f) {
        for (; smp>=1.0f; smp-=1.0f)
          cell = cell.walkDownDipFrom(y);
      } else {
        for (; smp<=-1.0f; smp+=1.0f)
          cell = cell.walkUpDipFrom(y);
      }

      // Account for any remaining fractional shift.
      float y1 = y[0], y2 = y[1], y3 = y[2];
      y1 += smp;
      y2 += smp*cell.us*cell.u2;
      y3 += smp*cell.us*cell.u3;

      // Unit dip vector.
      u1 = cell.u1;
      u2 = cell.u2;
      u3 = cell.u3;

      // Offset vector d.
      d1 = 0.0f;
      d2 =  _offset*cell.v3;
      d3 = -_offset*cell.v2;

      // Reflector slopes at point y+d.
      p2 = imageValueAt(y1+d1,y2+d2,y3+d3,_p2);
      p3 = imageValueAt(y1+d1,y2+d2,y3+d3,_p3);

      // Unit-vector a normal to reflector.
      a1 = 1.0f/sqrt(1.0f+p2*p2+p3*p3);
      a2 = -p2*a1;
      a3 = -p3*a1;

      // Slip adjustment for reflector slope on plus side.
      float ap = (a1*d1+a2*d2+a3*d3)/(a1*u1+a2*u2+a3*u3);
      float u1p = ap*u1;
      float u2p = ap*u2;
      float u3p = ap*u3;

      // Record total dip slip in cell at which we began the walk.
      cellBegin.s1 = y1-x1;
      cellBegin.s2 = y2-x2;
      cellBegin.s3 = y3-x3;
      if (!_zeroSlope) {
        cellBegin.s1 += u1m+u1p;
        cellBegin.s2 += u2m+u2p;
        cellBegin.s3 += u3m+u3p;
      }
    }
  }

  private static void clearErrors(FaultSkin skin) {
    for (FaultCell cell:skin)
      cell.emp = null;
  }

  private static float alignmentError(float f, float g) {
    float fmg = f-g;
    return fmg*fmg;
  }

  private static float imageValueAt(
    float p1, float p2, float p3, float[][][]f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    int i1 = max(0,min(n1-1,round(p1)));
    int i2 = max(0,min(n2-1,round(p2)));
    int i3 = max(0,min(n3-1,round(p3)));
    return f[i3][i2][i1];
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
