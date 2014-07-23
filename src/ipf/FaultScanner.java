/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package ipf;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Stopwatch;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.util.Parallel.*;

import static ipf.FaultGeometry.*;

/**
 * Computes fault likelihoods, strikes, and dips, by scanning over fault
 * orientations. Fault likelihoods are in the range [0,1], where 0 and 1
 * denote lowest and highest likelihoods, respectively. Computed fault strike
 * and dip angles are those for which maximum fault likelihoods occurred. 
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.06.29
 */
public class FaultScanner {

  /**
   * Constructs a fault scanner with specified parameters.
   * @param sigmaPhi half-width for smoothing along strike of fault planes.
   * @param sigmaTheta half-width for smoothing along dip of fault planes.
   */
  public FaultScanner(double sigmaPhi, double sigmaTheta) {
    _sigmaPhi = sigmaPhi;
    _sigmaTheta = sigmaTheta;
  }

  /**
   * Gets a sampling of fault strike phi appropriate for this scanner.
   * @param phiMin minimum fault strike, in degrees.
   * @param phiMax maximum fault strike, in degrees.
   */
  public Sampling getPhiSampling(double phiMin, double phiMax) {
    return angleSampling(_sigmaPhi,phiMin,phiMax);
  }

  /**
   * Gets a sampling of fault dip theta appropriate for this scanner.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   */
  public Sampling getThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
  }

  /**
   * Gets the frequencies of fault values (e.g., strikes or dips) in an array.
   * Each element of the returned array corresponds to a sampled value, and
   * contains the fraction of values in the specified array that are nearest
   * to that sampled value. In other words, the returned array is like a
   * histogram, but normalized so that the sum of all frequencies is one.
   * <p>
   * Fault likelihoods, if specified by a non-null array, may be used to
   * weight the counting of values in the array of values to be counted.
   * @param sv sampling of fault values in the returned array.
   * @param fv array of fault values to be counted.
   * @param fl array of fault likelihoods; null, for no weighting.
   * @return array of frequencies of fault values.
   */
  public static float[] getFrequencies(
    Sampling sv, float[][][] fv, float[][][] fl) {
    int n1 = fv[0][0].length;
    int n2 = fv[0].length;
    int n3 = fv.length;
    int nv = sv.getCount();
    float[] vf = new float[nv];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fvi = fv[i3][i2][i1];
          float fli = (fl!=null)?fl[i3][i2][i1]:1.0f;
          int iv = sv.indexOfNearest(fvi);
          if (iv>=0 && iv<nv)
            vf[iv] += fli;
        }
      }
    }
    float vfsum = sum(vf);
    float vfscl = (vfsum>0.0f)?1.0f/vfsum:1.0f;
    return mul(vf,vfscl);
  }

  /**
   * Returns slopes and planarities of features in a specified image.
   * Image features are assumed to be locally planar, with slopes that may
   * vary throughout the image. Smoothing parameters control the extents of
   * Gaussian windows within which slope is estimated for each image sample.
   * Typically, because slopes tend to vary most rapidly in the 2nd or 3rd
   * dimensions, the extent of smoothing for the 1st dimension should be
   * greater than that for the 2nd and 3rd dimensions.
   * @param sigma1 half-width for smoothing in 1st dimension.
   * @param sigma2 half-width for smoothing in 2nd dimension.
   * @param sigma3 half-width for smoothing in 3rd dimension.
   * @param slopeMax upper bound on computed slopes.
   * @param f input image for which to compute slopes.
   * @return array {p2,p3,ep} of slopes and planarities.
   */
  public static float[][][][] slopes(
      double sigma1, double sigma2, double sigma3, 
      double slopeMax, float[][][] f) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float p2min = (float)(-slopeMax);
    final float p2max = (float)( slopeMax);
    final float p3min = (float)(-slopeMax);
    final float p3max = (float)( slopeMax);

    // Normal vectors.
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    final float[][][] ep = new float[n3][n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(sigma1,sigma2,sigma3);
    lof.applyForNormalPlanar(f,u1,u2,u3,ep);

    // Slopes from normal vectors.
    final float[][][] p2 = u2;
    final float[][][] p3 = u3;
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
          if (-u2i<p2min*u1i) u2i = -p2min*u1i;
          if (-u2i>p2max*u1i) u2i = -p2max*u1i;
          if (-u3i<p3min*u1i) u3i = -p3min*u1i;
          if (-u3i>p3max*u1i) u3i = -p3max*u1i;
          if (u1i==0.0f) {
            p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
            p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
          } else {
            p2[i3][i2][i1] = -u2i/u1i;
            p3[i3][i2][i1] = -u3i/u1i;
          }
        }
      }
    }});
    return new float[][][][]{p2,p3,ep};
  }

  /**
   * Returns an image with specified samples taper to zero at edges.
   * Tapering enables simple zero-value boundary conditions to be
   * used in fault scanning without artifacts caused by abrupt 
   * truncations of strong image features near edges.
   * @param m1 width of the tapered samples near edges in 1st dimension.
   * @param m2 width of the tapered samples near edges in 2nd dimension.
   * @param m3 width of the tapered samples near edges in 3rd dimension.
   * @param f input image.
   * @return the tapered image.
   */
  public static float[][][] taper(int m1, int m2, int m3, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] g = copy(f);
    float[] t1 = new float[m1];
    float[] t2 = new float[m2];
    float[] t3 = new float[m3];
    for (int i1=0; i1<m1; ++i1)
      t1[i1] = (float)(0.54+0.46*cos(PI*(m1-i1)/m1));
    for (int i2=0; i2<m2; ++i2)
      t2[i2] = (float)(0.54+0.46*cos(PI*(m2-i2)/m2));
    for (int i3=0; i3<m3; ++i3)
      t3[i3] = (float)(0.54+0.46*cos(PI*(m3-i3)/m3));
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0,j1=n1-1; i1<m1; ++i1,--j1) {
          float ti = t1[i1];
          g[i3][i2][i1] *= ti;
          g[i3][i2][j1] *= ti;
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0,j2=n2-1; i2<m2; ++i2,--j2) {
        float ti = t2[i2];
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[i3][j2][i1] *= ti;
        }
      }
    }
    for (int i3=0,j3=n3-1; i3<m3; ++i3,--j3) {
      float ti = t3[i3];
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          g[i3][i2][i1] *= ti;
          g[j3][i2][i1] *= ti;
        }
      }
    }
    return g;
  }

  /**
   * Scans a specified image for fault strikes and dips.
   * @param phiMin minimum fault strike, in degrees.
   * @param phiMax maximum fault strike, in degrees.
   * @param thetaMin minimum fault dip, in degrees.
   * @param thetaMax maximum fault dip, in degrees.
   * @param p2 slopes in the 2nd dimension.
   * @param p3 slopes in the 3rd dimension.
   * @param g the image to be scanned.
   * @return array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   */
  public float[][][][] scan(
      double phiMin, double phiMax,
      double thetaMin, double thetaMax,
      float[][][] p2, float[][][] p3, float[][][] g) {
    Sampling sp = makePhiSampling(phiMin,phiMax);
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    return scan(sp,st,p2,p3,g);
  }

  /**
   * Scans with the specified sampling of fault strikes and dips.
   * @param phiSampling sampling of fault strikes, in degrees.
   * @param thetaSampling sampling of fault dip angles, in degrees.
   * @param p2 slopes in the 2nd dimension.
   * @param p3 slopes in the 3rd dimension.
   * @param g the image to be scanned.
   * @return array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   */
  public float[][][][] scan(
      Sampling phiSampling, Sampling thetaSampling,
      float[][][] p2, float[][][] p3, float[][][] g) {
    float[][][][] snd = semblanceNumDen(p2,p3,g);
    return scan(phiSampling,thetaSampling,snd);
  }

  /**
   * Thins fault images to include only ridges in fault likelihoods.
   * After thinning, may be only one voxel wide. Thinned fault strikes and
   * dips are set to zero where thinned fault likelihoods are zero.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   * @return array {fl,fp,ft} of thinned fault likelihoods, strikes, and dips.
   */
  public static float[][][][] thin(float[][][][] flpt) {
    int n1 = n1(flpt);
    int n2 = n2(flpt);
    int n3 = n3(flpt);
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    f = copy(f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.applyX0X(f,f);
    rgf.applyXX0(f,f);
    float[][][] ff = new float[n3][n2][n1];
    float[][][] pp = new float[n3][n2][n1];
    float[][][] tt = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p00 = p[i3 ][i2 ];
        float[] t00 = t[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          float f000 = f00[i1];
          float p000 = p00[i1];
          float t000 = t00[i1];
          if ((                p000<= 22.5f && f0m[i1]<f000 && f0p[i1]<f000) ||
              ( 22.5f<=p000 && p000<= 67.5f && fpm[i1]<f000 && fmp[i1]<f000) ||
              ( 67.5f<=p000 && p000<=112.5f && fp0[i1]<f000 && fm0[i1]<f000) ||
              (112.5f<=p000 && p000<=157.5f && fpp[i1]<f000 && fmm[i1]<f000) ||
              (157.5f<=p000 && p000<=202.5f && f0p[i1]<f000 && f0m[i1]<f000) ||
              (202.5f<=p000 && p000<=247.5f && fmp[i1]<f000 && fpm[i1]<f000) ||
              (247.5f<=p000 && p000<=292.5f && fm0[i1]<f000 && fp0[i1]<f000) ||
              (292.5f<=p000 && p000<=337.5f && fmm[i1]<f000 && fpp[i1]<f000) ||
              (337.5f<=p000                 && f0m[i1]<f000 && f0p[i1]<f000)) {
            ff[i3][i2][i1] = f000;
            pp[i3][i2][i1] = p000;
            tt[i3][i2][i1] = t000;
          } else {
            pp[i3][i2][i1] = NO_STRIKE;
            tt[i3][i2][i1] = NO_DIP;
          }
        }
      }
    }
    float[][][][] flptn = new float[][][][]{ff,pp,tt};
    removeEdgeEffects(flptn);
    return flptn;
  }

  /**
   * Applies structure-oriented smoothing limited by fault likelihoods.
   * For this method, faults are assumed to exist and smoothing stops
   * at samples where fault likelihood exceeds a specified value.
   * This method is usually applied using thinned fault likelihoods.
   * @param flstop smoothing stops where fault likelihood &gt; this value.
   * @param sigma smoothing radius (except near faults).
   * @param p2 array of slopes in 2nd dimension.
   * @param p3 array of slopes in 3rd dimension.
   * @param fl array of fault likelihoods, typically thinned.
   * @param g image to be smoothed.
   */
  public static float[][][] smooth(
      double flstop, double sigma, float[][][] p2, float[][][] p3, 
      float[][][] fl, float[][][] g) {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    EigenTensors3 d = new EigenTensors3(n1,n2,n3,true);
    d.setEigenvalues(0.001f,1.00f,1.00f);
    float[][][] s = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          s[i3][i2][i1] = (fl[i3][i2][i1]<flstop)?1.0f:0.0f;
          float p2i = p2[i3][i2][i1];
          float p3i = p3[i3][i2][i1];
          float u1i = 1.0f/sqrt(1.0f+p2i*p2i+p3i*p3i);
          float u2i = -p2i*u1i;
          float u3i = -p3i*u1i;
          float usi = 1.0f/sqrt(u1i*u1i+u2i*u2i);
          float w1i = -u2i*usi;
          float w2i =  u1i*usi;
          float w3i = 0.0f;
          d.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
          d.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
        }
      }
    }
    float c = (float)(0.5*sigma*sigma);
    float[][][] h = new float[n3][n2][n1];
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(d,c,s,g,h);
    return h;
  }

  /**
   * Adjusts dips for a specified aspect ratio dz/dx. A fault scanner computes
   * fault dips in degrees measured in sample coordinates. Because vertical
   * image sampling intervals are often less than horizontal sampling
   * intervals, fault dips measured in sample coordinates tend to be greater
   * than those in physical coordinates. This method converts fault dips to
   * degrees measured in physical coordinates, using the specified ratio of
   * vertical-to-horizontal sampling intervals.
   * @param dzdx ratio of vertical to horizontal sampling intervals.
   * @param ft array of fault dips measured in sample coordinates.
   * @return array of fault dips measured in physical coordinates.
   */
  public static float[][][] convertDips(double dzdx, float[][][] ft) {
    float scale = (float)dzdx;
    int n1 = ft[0][0].length;
    int n2 = ft[0].length;
    int n3 = ft.length;
    float[][][] gt = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fti = ft[i3][i2][i1];
          if (fti!=NO_DIP)
            gt[i3][i2][i1] = toDegrees(atan(scale*tan(toRadians(fti))));
        }
      }
    }
    return gt;
  }

  /**
   * Adjust strikes for coordinate system and rotation.
   * @param lh true, if left-handed coordinate system; false, otherwise.
   * @param pa, amount to add to azimuths.
   * @param fp array of fault strikes measured in sample coordinates.
   * @return array of fault strikes measured in physical coordinates.
   */
  public static float[][][] convertStrikes(
      boolean lh, double pa, float[][][] fp) {
    int n1 = fp[0][0].length;
    int n2 = fp[0].length;
    int n3 = fp.length;
    float fpa = (float)pa;
    float[][][] gp = copy(fp);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fpi = fp[i3][i2][i1];
          if (fpi!=NO_STRIKE)
            gp[i3][i2][i1] = range360((lh?360.0f-fpi:fpi)+fpa);
        }
      }
    }
    return gp;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sigmaPhi,_sigmaTheta;

  private static final float NO_STRIKE = -0.00001f;
  private static final float NO_DIP    = -0.00001f;

  private static void trace(String s) {
    System.out.println(s);
  }

  // This scan smooths semblance numerators and denominators along fault
  // planes by first rotating and shearing those images before applying
  // fast recursive axis-aligned smoothing filters.
  private float[][][][] scan(
      Sampling phiSampling, Sampling thetaSampling,
      float[][][][] snd) {
    // Algorithm: given snum,sden (semblance numerators and denominators)
    // initialize f,p,t (fault likelihood, phi, and theta)
    // for all phi:
    //   rotate snum,sden so that strike vector is aligned with axis 2
    //   smooth snum,sden along fault strike (that is, along axis 2)
    //   compute fphi,tphi (fault likelihood and dip) in 1-3 slices
    //   unrotate fphi,tphi to original coordinates
    //   update f,p,t for maximum likelihood
    final int n1 = snd[0][0][0].length;
    final int n2 = snd[0][0].length;
    final int n3 = snd[0].length;
    final float[][][] f = new float[n3][n2][n1];
    final float[][][] p = new float[n3][n2][n1];
    final float[][][] t = new float[n3][n2][n1];
    final float tmin = (float)thetaSampling.getFirst();
    final float tmax = (float)thetaSampling.getLast();
    int np = phiSampling.getCount();
    Stopwatch sw = new Stopwatch();
    sw.start();
    for (int ip=0; ip<np; ++ip) {
      final float phi = (float)phiSampling.getValue(ip);
      if (ip>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)np/(double)ip-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("FaultScanner.scan: done in "+timeLeftSec+" seconds");
      }
      Rotator r = new Rotator(phi,n1,n2,n3);
      float[][][][] rsnd = r.rotate(snd);
      smooth2(rsnd);
      float[][][][] rftp = scanTheta(thetaSampling,rsnd);
      rsnd = null; // enable gc to collect this large array
      float[][][][] ftp = r.unrotate(rftp);
      final float[][][] fp = ftp[0];
      final float[][][] tp = ftp[1];
      loop(n3,new LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] f32 = f[i3][i2];
          float[] p32 = p[i3][i2];
          float[] t32 = t[i3][i2];
          float[] fp32 = fp[i3][i2];
          float[] tp32 = tp[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float fpi = fp32[i1];
            float tpi = tp32[i1];
            if (fpi<0.0f) fpi = 0.0f; // necessary because of sinc
            if (fpi>1.0f) fpi = 1.0f; // interpolation in unrotate,
            if (tpi<tmin) tpi = tmin; // for both fault likelihood
            if (tpi>tmax) tpi = tmax; // and fault dip theta
            if (fpi>f32[i1]) {
              f32[i1] = fpi;
              p32[i1] = phi;
              t32[i1] = tpi;
            }
          }
        }
      }});
    }
    sw.stop();
    trace("FaultScanner.scan: done");
    return new float[][][][]{f,p,t};
  }

  // Sampling of angles depends on extent of smoothing.
  private Sampling makePhiSampling(double phiMin, double phiMax) {
    return angleSampling(_sigmaPhi,phiMin,phiMax);
  }
  private Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigmaTheta,thetaMin,thetaMax);
  }
  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  // Numbers of samples in 3D arrays (arrays of arrays of arrays),
  // which after rotation may contain some null arrays.
  private static int n1(float[][][] f) {
    int n1 = 0;
    int n2 = f[0].length;
    int n3 = f.length;
    for (int i3=0; i3<n3 && n1==0; ++i3) {
      for (int i2=0; i2<n2 && n1==0; ++i2) {
        if (f[i3][i2]!=null)
          n1 = f[i3][i2].length;
      }
    }
    return n1;
  }
  private static int n2(float[][][] f) {
    return f[0].length;
  }
  private static int n3(float[][][] f) {
    return f.length;
  }
  private static int n1(float[][][][] f) {
    return n1(f[0]);
  }
  private static int n2(float[][][][] f) {
    return n2(f[0]);
  }
  private static int n3(float[][][][] f) {
    return n3(f[0]);
  }

  // Get/set non-null slices of rotated 3D arrays
  private static float[][] extractSlice2(int i2, float[][][] x) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    float[][] x2 = (m3>0)?new float[m3][n1]:null;
    for (int i3=0; i3<m3; ++i3)
      copy(x[i3+i3lo][i2],x2[i3]);
    return x2;
  }
  private static void restoreSlice2(int i2, float[][][] x, float[][] x2) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i3lo = i3lo(i2,x);
    int i3hi = i3hi(i2,x);
    int m3 = 1+i3hi-i3lo;
    assert x2.length==m3:"x2 length is correct";
    for (int i3=0; i3<m3; ++i3)
      copy(x2[i3],x[i3+i3lo][i2]);
  }
  private static float[][] extractSlice3(int i3, float[][][] x) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i2lo = i2lo(i3,x);
    int i2hi = i2hi(i3,x);
    int m2 = 1+i2hi-i2lo;
    float[][] x3 = (m2>0)?new float[m2][n1]:null;
    for (int i2=0; i2<m2; ++i2)
      copy(x[i3][i2+i2lo],x3[i2]);
    return x3;
  }
  private static void restoreSlice3(int i3, float[][][] x, float[][] x3) {
    int n1 = n1(x);
    int n2 = n2(x);
    int n3 = n3(x);
    int i2lo = i2lo(i3,x);
    int i2hi = i2hi(i3,x);
    int m2 = 1+i2hi-i2lo;
    assert x3.length==m2:"x3 length is correct";
    for (int i2=0; i2<m2; ++i2)
      copy(x3[i2],x[i3][i2+i2lo]);
  }
  private static int i2lo(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2lo = 0;
    while (i2lo<n2 && x[i3][i2lo]==null)
      ++i2lo;
    return i2lo;
  }
  private static int i2hi(int i3, float[][][] x) {
    int n2 = x[0].length;
    int i2hi = n2-1;
    while (i2hi>=0 && x[i3][i2hi]==null)
      --i2hi;
    return i2hi;
  }
  private static int i3lo(int i2, float[][][] x) {
    int n3 = x.length;
    int i3lo = 0;
    while (i3lo<n3 && x[i3lo][i2]==null)
      ++i3lo;
    return i3lo;
  }
  private static int i3hi(int i2, float[][][] x) {
    int n3 = x.length;
    int i3hi = n3-1;
    while (i3hi>=0 && x[i3hi][i2]==null)
      --i3hi;
    return i3hi;
  }

  // Shear horizontally such that q(i1,i2) = p(i1,i2+s*i1).
  private static float[][] shear(
    SincInterpolator si, double s, float[][] p)
  {
    int n1 = p[0].length;
    int n2p = p.length;
    int n2q = n2p+(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] q = new float[n2q][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2p; ++i2)
        pp[i2] = p[i2][i1];
      double f2q = (s<0.0f)?s*i1:s*i1-dqp;
      si.interpolate(n2p,1.0,0.0,pp,n2q,1.0f,f2q,qq);
      for (int i2=0; i2<n2q; ++i2)
        q[i2][i1] = qq[i2];
    }
    return q;
  }

  // Unshear horizontally such that p(i1,i2) = q(i1,i2-s*i1).
  private static float[][] unshear(
    SincInterpolator si, double s, float[][] q)
  {
    int n1 = q[0].length;
    int n2q = q.length;
    int n2p = n2q-(int)(abs(s)*n1);
    double dqp = n2q-n2p;
    float[][] p = new float[n2p][n1]; 
    float[] pp = new float[n2p];
    float[] qq = new float[n2q];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2q; ++i2)
        qq[i2] = q[i2][i1];
      double f2p = (s<0.0f)?-s*i1:-s*i1+dqp;
      si.interpolate(n2q,1.0,0.0,qq,n2p,1.0f,f2p,pp);
      for (int i2=0; i2<n2p; ++i2)
        p[i2][i1] = pp[i2];
    }
    return p;
  }

  // Horizontal smoothing of rotated snum,sden along axis 2.
  private void smooth2(final float[][][][] snd) {
    final int n1 = n1(snd), n2 = n2(snd), n3 = n3(snd);
    final RecursiveExponentialFilter ref = makeRef(_sigmaPhi);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      for (int is=0; is<2; ++is) {
        float[][] s3 = extractSlice3(i3,snd[is]);
        if (s3!=null) {
          ref.apply2(s3,s3); 
          restoreSlice3(i3,snd[is],s3);
        }
      }
    }});
  }

  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // For one fault strike, scans over all fault dips theta. The fault strike
  // vector has already been aligned with image axis 2, and semblance
  // numerators and denominators have already been smoothed in that direction.
  // Therefore, this scan over fault dip can be performed independently for
  // each i2. For each fault dip theta, this method shears semblance num and
  // den to align any faults having that dip with image axis 1. Semblance
  // num and den are then smoothed vertically, with an extent sigma that is
  // dip-adjusted (shorter for smaller fault dips), so that after unshearing
  // the extent of smoothing is roughly the same for all fault dips.
  private float[][][][] scanTheta(Sampling thetaSampling, float[][][][] snd) {
    final int n1 = n1(snd), n2 = n2(snd), n3 = n3(snd);
    final Sampling st = thetaSampling;
    final float[][][][] ft = like(snd);
    final float[][][] sn = snd[0];
    final float[][][] sd = snd[1];
    final float[][][] f = ft[0];
    final float[][][] t = ft[1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n2,new LoopInt() {
    public void compute(int i2) {
      float[][] sn2 = extractSlice2(i2,sn);
      float[][] sd2 = extractSlice2(i2,sd);
      if (sn2==null)
        return;
      int n3 = sn2.length;
      int nt = st.getCount();
      for (int it=0; it<nt; ++it) {
        float ti = (float)st.getValue(it);
        float theta = toRadians(ti);
        float shear = -1.0f/tan(theta);
        float[][] sns = shear(si,shear,sn2);
        float[][] sds = shear(si,shear,sd2);
        float sigma = (float)_sigmaTheta*sin(theta);
        RecursiveExponentialFilter ref = makeRef(sigma);
        ref.apply1(sns,sns);
        ref.apply1(sds,sds);
        float[][] ss = semblanceFromNumDen(sns,sds);
        float[][] s2 = unshear(si,shear,ss);
        for (int i3=0,j3=i3lo(i2,f); i3<n3; ++i3,++j3) {
          float[] s32 = s2[i3];
          float[] f32 = f[j3][i2];
          float[] t32 = t[j3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float st = s32[i1]; // semblance
            st = st*st; // semblance^2
            st = st*st; // semblance^4
            st = st*st; // semblance^8
            float fi = 1.0f-st;
            if (fi>f32[i1]) {
              f32[i1] = fi;
              t32[i1] = ti;
            }
          }
        }
      }
    }});
    return ft;
  }

  // Makes an array like that specified, including any null arrays.
  private float[][][][] like(float[][][][] p) {
    int n1 = n1(p[0]);
    int n2 = n2(p[0]);
    int n3 = n3(p[0]);
    int np = p.length;
    float[][][][] q = new float[np][n3][n2][];
    for (int ip=0; ip<np; ++ip) {
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          q[ip][i3][i2] = (p[ip][i3][i2]!=null)?new float[n1]:null;
        }
      }
    }
    return q;
  }

  // Computes fault semblance numerators and denominators.
  private static float[][][][] semblanceNumDen(
    final float[][][] p2, final float[][][] p3, final float[][][] f) 
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] sn = new float[n3][n2][n1];
    final float[][][] sd = new float[n3][n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    loop(n3,new LoopInt() {
    public void compute(int i3) {
      float[] xmm = new float[n1];
      float[] xm0 = new float[n1];
      float[] xmp = new float[n1];
      float[] x0m = new float[n1];
      float[] x0p = new float[n1];
      float[] xpm = new float[n1];
      float[] xp0 = new float[n1];
      float[] xpp = new float[n1];
      float[] gmm = new float[n1];
      float[] gm0 = new float[n1];
      float[] gmp = new float[n1];
      float[] g0m = new float[n1];
      float[] g0p = new float[n1];
      float[] gpm = new float[n1];
      float[] gp0 = new float[n1];
      float[] gpp = new float[n1];
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p2mm = p2[i3m][i2m];
        float[] p2m0 = p2[i3m][i2 ];
        float[] p2mp = p2[i3m][i2p];
        float[] p20m = p2[i3 ][i2m];
        float[] p200 = p2[i3 ][i2 ];
        float[] p20p = p2[i3 ][i2p];
        float[] p2pm = p2[i3p][i2m];
        float[] p2p0 = p2[i3p][i2 ];
        float[] p2pp = p2[i3p][i2p];
        float[] p3mm = p3[i3m][i2m];
        float[] p3m0 = p3[i3m][i2 ];
        float[] p3mp = p3[i3m][i2p];
        float[] p30m = p3[i3 ][i2m];
        float[] p300 = p3[i3 ][i2 ];
        float[] p30p = p3[i3 ][i2p];
        float[] p3pm = p3[i3p][i2m];
        float[] p3p0 = p3[i3p][i2 ];
        float[] p3pp = p3[i3p][i2p];
        float[] sn32 = sn[i3][i2];
        float[] sd32 = sd[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xmm[i1] = i1-p3mm[i1]-p2mm[i1];
          xm0[i1] = i1-p3m0[i1]         ;
          xmp[i1] = i1-p3mp[i1]+p2mp[i1];
          x0m[i1] = i1         -p20m[i1];
          x0p[i1] = i1         +p20p[i1];
          xpm[i1] = i1+p3pm[i1]-p2pm[i1];
          xp0[i1] = i1+p3p0[i1]         ;
          xpp[i1] = i1+p3pp[i1]+p2pp[i1];
        }
        si.interpolate(n1,1.0,0.0,fmm,n1,xmm,gmm);
        si.interpolate(n1,1.0,0.0,fm0,n1,xm0,gm0);
        si.interpolate(n1,1.0,0.0,fmp,n1,xmp,gmp);
        si.interpolate(n1,1.0,0.0,f0m,n1,x0m,g0m);
        si.interpolate(n1,1.0,0.0,f0p,n1,x0p,g0p);
        si.interpolate(n1,1.0,0.0,fpm,n1,xpm,gpm);
        si.interpolate(n1,1.0,0.0,fp0,n1,xp0,gp0);
        si.interpolate(n1,1.0,0.0,fpp,n1,xpp,gpp);
        float[] hmm = gmm, hm0 = gm0, hmp = gmp;
        float[] h0m = g0m, h00 = f00, h0p = g0p;
        float[] hpm = gpm, hp0 = gp0, hpp = gpp;
        if (            i2==0   ) h0m = h00;
        if (            i2==n2-1) h0p = h00;
        if (i3==0               ) hm0 = h00;
        if (i3==n3-1            ) hp0 = h00;
        if (i3==0    && i2==0   ) hmm = h00;
        if (i3==0    && i2==n2-1) hmp = h00;
        if (i3==n3-1 && i2==0   ) hpm = h00;
        if (i3==n3-1 && i2==n2-1) hpp = h00;
        for (int i1=0; i1<n1; ++i1) {
          float hmmi = hmm[i1];
          float hm0i = hm0[i1];
          float hmpi = hmp[i1];
          float h0mi = h0m[i1];
          float h00i = h00[i1];
          float h0pi = h0p[i1];
          float hpmi = hpm[i1];
          float hp0i = hp0[i1];
          float hppi = hpp[i1];
          float sumn = hmmi+hm0i+hmpi+
                       h0mi+h00i+h0pi+
                       hpmi+hp0i+hppi;
          float sumd = hmmi*hmmi+hm0i*hm0i+hmpi*hmpi+
                       h0mi*h0mi+h00i*h00i+h0pi*h0pi+
                       hpmi*hpmi+hp0i*hp0i+hppi*hppi;
          sn32[i1] = sumn*sumn;
          sd32[i1] = 9.0f*sumd;
        }
      }
    }});
    return new float[][][][]{sn,sd};
  }

  // Computes semblance ratios from numerators and denominators.
  // Takes care to ensure that semblances are in range [0,1].
  private static float[][] semblanceFromNumDen(float[][] sn, float[][] sd) {
    int n1 = sn[0].length;
    int n2 = sn.length;
    float[][] sr = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] sn2 = sn[i2];
      float[] sd2 = sd[i2];
      float[] sr2 = sr[i2];
      for (int i1=0; i1<n1; ++i1) {
        float sni = sn2[i1];
        float sdi = sd2[i1];
        if (sdi<=0.0f || sni<=0.0f) {
          sr2[i1] = 0.0f;
        } else if (sdi<sni) {
          sr2[i1] = 1.0f;
        } else {
          sr2[i1] = sni/sdi;
        }
      }
    }
    return sr;
  }

  // Removes spurious faults caused by image boundaries. A sample of fault
  // likelihood, strike and dip is deemed spurious if it is both near and
  // nearly parallel to the image boundary. This method zeros likelihoods,
  // strikes and dips for any such samples.
  private static void removeEdgeEffects(float[][][][] flpt) {
    int n1 = n1(flpt);
    int n2 = n2(flpt);
    int n3 = n3(flpt);
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int imax = 5; // max number of samples near boundary
    float amin = 30.0f; // min angle between normal vectors
    float cmax = cos(toRadians(amin));
    float wwmax = cmax*cmax; 
    for (int i3=0,j3=n3-1; i3<imax; ++i3,--j3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fi = f[i3][i2][i1];
          float pi = p[i3][i2][i1];
          float ti = t[i3][i2][i1];
          float fj = f[j3][i2][i1];
          float pj = p[j3][i2][i1];
          float tj = t[j3][i2][i1];
          if (fi!=0.0f) {
            float[] w = faultNormalVectorFromStrikeAndDip(pi,ti);
            float w3 = w[2];
            if (w3*w3>wwmax) {
              f[i3][i2][i1] = 0.0f;
              p[i3][i2][i1] = 0.0f;
              t[i3][i2][i1] = 0.0f;
            }
          }
          if (fj!=0.0f) {
            float[] w = faultNormalVectorFromStrikeAndDip(pj,tj);
            float w3 = w[2];
            if (w3*w3>wwmax) {
              f[j3][i2][i1] = 0.0f;
              p[j3][i2][i1] = 0.0f;
              t[j3][i2][i1] = 0.0f;
            }
          }
        }
      }
    }
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0,j2=n2-1; i2<imax; ++i2,--j2) {
        for (int i1=0; i1<n1; ++i1) {
          float fi = f[i3][i2][i1];
          float pi = p[i3][i2][i1];
          float ti = t[i3][i2][i1];
          float fj = f[i3][j2][i1];
          float pj = p[i3][j2][i1];
          float tj = t[i3][j2][i1];
          if (fi!=0.0f) {
            float[] w = faultNormalVectorFromStrikeAndDip(pi,ti);
            float w2 = w[1];
            if (w2*w2>wwmax) {
              f[i3][i2][i1] = 0.0f;
              p[i3][i2][i1] = 0.0f;
              t[i3][i2][i1] = 0.0f;
            }
          }
          if (fj!=0.0f) {
            float[] w = faultNormalVectorFromStrikeAndDip(pj,tj);
            float w2 = w[1];
            if (w2*w2>wwmax) {
              f[i3][j2][i1] = 0.0f;
              p[i3][j2][i1] = 0.0f;
              t[i3][j2][i1] = 0.0f;
            }
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // image rotator

  private static class Rotator {

    Rotator(double phi, int n1, int n2, int n3) {
      _n1 = n1;

      // angle phi in radians, cosine and sine
      _phir = toRadians(phi);
      _cosp = cos(_phir);
      _sinp = sin(_phir);

      // center of rotation
      _x2c = 0.5*(n2-1.0);
      _x3c = 0.5*(n3-1.0);

      // input sampling
      _s2p = new Sampling(n2,1.0,0.0);
      _s3p = new Sampling(n3,1.0,0.0);

      // corners of input sampling rectangle
      double[] x2s = { 0.0, 0.0,n2-1,n2-1};
      double[] x3s = { 0.0,n3-1,n3-1, 0.0};

      // bounds after rotation
      double x2min =  Double.MAX_VALUE;
      double x3min =  Double.MAX_VALUE;
      double x2max = -Double.MAX_VALUE;
      double x3max = -Double.MAX_VALUE;
      for (int i=0; i<4; ++i) {
        double x2q = x2q(x2s[i],x3s[i]);
        double x3q = x3q(x2s[i],x3s[i]);
        if (x2q<x2min) x2min = x2q;
        if (x2q>x2max) x2max = x2q;
        if (x3q<x3min) x3min = x3q;
        if (x3q>x3max) x3max = x3q;
      }
      x2min = floor(x2min);
      x2max = ceil(x2max);
      x3min = floor(x3min);
      x3max = ceil(x3max);

      // sampling after rotation
      int n2q = max(2,1+(int)(x2max-x2min+0.5));
      int n3q = max(2,1+(int)(x3max-x3min+0.5));
      double d2q = 1.0;
      double d3q = 1.0;
      double f2q = x2min;
      double f3q = x3min;
      _s2q = new Sampling(n2q,d2q,f2q);
      _s3q = new Sampling(n3q,d3q,f3q);
      //trace("s2p: n2p="+n2);
      //trace("s3p: n3p="+n3);
      //trace("s2q: n2q="+n2q+" d2q="+d2q+" f2q="+f2q);
      //trace("s3q: n3q="+n3q+" d3q="+d3q+" f3q="+f3q);
    }

    float[][][][] rotate(float[][][][] p) {
      int n = p.length;
      float[][][][] q = new float[n][][][];
      for (int i=0; i<n; ++i)
        q[i] = rotate(p[i]);
      return q;
    }

    float[][][][] unrotate(float[][][][] p) {
      int n = p.length;
      float[][][][] q = new float[n][][][];
      for (int i=0; i<n; ++i)
        q[i] = unrotate(p[i]);
      return q;
    }

    float[][][] rotate(float[][][] p) {
      final float[][][] fp = p;
      final float[][] siTable = _siTable;
      final int nsinc = siTable.length;
      final int lsinc = siTable[0].length;
      final Sampling s2p = _s2p;
      final Sampling s3p = _s3p;
      final Sampling s2q = _s2q;
      final Sampling s3q = _s3q;
      final int n1 = _n1;
      final int n2p = _s2p.getCount();
      final int n3p = _s3p.getCount();
      final int n2q = _s2q.getCount();
      final int n3q = _s3q.getCount();
      final float[][][] q = new float[n3q][n2q][];
      loop(n3q,new LoopInt() {
        public void compute(int i3) {
          double x3q = s3q.getValue(i3);
          for (int i2=0; i2<n2q; ++i2) {
            double x2q = s2q.getValue(i2);
            double x2p = x2p(x2q,x3q);
            double x3p = x3p(x2q,x3q);
            if (inBounds(x2p,x3p)) {
              float[] q32 = q[i3][i2] = new float[n1];
              int i2p = (int)floor(x2p);
              int i3p = (int)floor(x3p);
              double f2p = x2p-i2p;
              double f3p = x3p-i3p;
              int k2p = (int)(f2p*(nsinc-1)+0.5);
              int k3p = (int)(f3p*(nsinc-1)+0.5);
              for (int k3s=0; k3s<lsinc; ++k3s) {
                float s3 = siTable[k3p][k3s];
                int j3p = i3p+k3s-lsinc/2+1;
                if (j3p<   0) j3p = 0;
                if (j3p>=n3p) j3p = n3p-1;
                for (int k2s=0; k2s<lsinc; ++k2s) {
                  float s2 = siTable[k2p][k2s];
                  int j2p = i2p+k2s-lsinc/2+1;
                  if (j2p<   0) j2p = 0;
                  if (j2p>=n2p) j2p = n2p-1;
                  float[] p32 = fp[j3p][j2p];
                  float s32 = s3*s2;
                  for (int i1=0; i1<n1; ++i1)
                    q32[i1] += p32[i1]*s32;
                }
              }
            }
          }
        }
      });
      return q;
    }

    float[][][] unrotate(float[][][] q) {
      final float[][][] fq = q;
      final float[][] siTable = _siTable;
      final int nsinc = siTable.length;
      final int lsinc = siTable[0].length;
      final Sampling s2p = _s2p;
      final Sampling s3p = _s3p;
      final Sampling s2q = _s2q;
      final Sampling s3q = _s3q;
      final int n1 = _n1;
      final int n2p = s2p.getCount();
      final int n3p = s3p.getCount();
      final int n2q = s2q.getCount();
      final int n3q = s3q.getCount();
      //System.out.println("n2p="+n2p+" n3p="+n3p+" n2q="+n2q+" n3q="+n3q);
      final double d2q = s2q.getDelta();
      final double d3q = s3q.getDelta();
      final double f2q = s2q.getFirst();
      final double f3q = s3q.getFirst();
      final float[][][] p = new float[n3p][n2p][n1];
      loop(n3p,new LoopInt() {
        public void compute(int i3) {
          double x3p = s3p.getValue(i3);
          for (int i2=0; i2<n2p; ++i2) {
            float[] p32 = p[i3][i2];
            double x2p = s2p.getValue(i2);
            double x2q = x2q(x2p,x3p);
            double x3q = x3q(x2p,x3p);
            double y2q = (x2q-f2q)/d2q;
            double y3q = (x3q-f3q)/d3q;
            int i2q = (int)floor(y2q);
            int i3q = (int)floor(y3q);
            double e2q = y2q-i2q;
            double e3q = y3q-i3q;
            int k2q = (int)(e2q*(nsinc-1)+0.5);
            int k3q = (int)(e3q*(nsinc-1)+0.5);
            for (int k3s=0; k3s<lsinc; ++k3s) {
              float s3 = siTable[k3q][k3s];
              int j3q = i3q+k3s-lsinc/2+1;
              if (j3q<   0) j3q = 0;
              if (j3q>=n3q) j3q = n3q-1;
              for (int k2s=0; k2s<lsinc; ++k2s) {
                float s2 = siTable[k2q][k2s];
                int j2q = i2q+k2s-lsinc/2+1;
                if (j2q<   0) j2q = 0;
                if (j2q>=n2q) j2q = n2q-1;
                float[] q32 = fq[j3q][j2q];
                if (q32!=null) {
                  float s32 = s3*s2;
                  for (int i1=0; i1<n1; ++i1)
                    p32[i1] += q32[i1]*s32;
                }
              }
            }
          }
        }
      });
      return p;
    }

    /////////////////////////////////////////////////////////////////////////
    // private

    private int _n1; // number of samples in 1st dimension
    private double _phir,_cosp,_sinp; // angle phi in radians, cosine, sine
    private double _x2c,_x3c; // coordinates of center of rotation
    private Sampling _s2p,_s3p; // samplings in original coordinates
    private Sampling _s2q,_s3q; // samplings in rotated coordinates
    private static float[][] _siTable; // sinc interpolation coefficients
    private static int HALF_LSINC; // half length of sinc interpolator
    static {
      SincInterpolator si = new SincInterpolator();
      _siTable = si.getTable();
      HALF_LSINC = _siTable[0].length/2;
    }
    private double x2p(double x2q, double x3q) {
      return _x2c+(x2q-_x2c)*_sinp-(x3q-_x3c)*_cosp;
    }
    private double x3p(double x2q, double x3q) {
      return _x3c+(x2q-_x2c)*_cosp+(x3q-_x3c)*_sinp;
    }
    private double x2q(double x2p, double x3p) {
      return _x2c+(x2p-_x2c)*_sinp+(x3p-_x3c)*_cosp;
    }
    private double x3q(double x2p, double x3p) {
      return _x3c-(x2p-_x2c)*_cosp+(x3p-_x3c)*_sinp;
    }
    /*
    private double x2p(double x2q, double x3q) {
      return _x2c+(x2q-_x2c)*_cosp-(x3q-_x3c)*_sinp;
    }
    private double x3p(double x2q, double x3q) {
      return _x3c+(x2q-_x2c)*_sinp+(x3q-_x3c)*_cosp;
    }
    private double x2q(double x2p, double x3p) {
      return _x2c+(x2p-_x2c)*_cosp+(x3p-_x3c)*_sinp;
    }
    private double x3q(double x2p, double x3p) {
      return _x3c-(x2p-_x2c)*_sinp+(x3p-_x3c)*_cosp;
    }
    */
    private boolean inBounds(double x2p, double x3p) {
      return _s2p.getFirst()-HALF_LSINC<=x2p && 
             _s3p.getFirst()-HALF_LSINC<=x3p && 
              x2p<=_s2p.getLast()+HALF_LSINC &&
              x3p<=_s3p.getLast()+HALF_LSINC;
    }
  }
}
