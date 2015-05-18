/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import java.util.logging.Logger;

/**
 * Iterative conjugate-gradient solver for Ax = b, where A is a square matrix.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class CgSolver {

  /**
   * The solver stops when one of these conditions is met.
   */
  public enum Stop {
    /**
     * Norm of residuals ||r|| = ||b-Ax|| is less than tiny * ||b||.
     * When using a preconditioner M, ||r|| = sqrt((b-Ax)'M(b-Ax)).
     */
    TINY,
    /**
     * The maximum number of iterations was performed.
     */
    MAXI,
  }

  /**
   * Information returned by this iterative solver.
   */
  public static class Info {
    private Info(Stop stop, int niter, double bnorm, double rnorm) {
      this.stop = stop;
      this.niter = niter;
      this.bnorm = bnorm;
      this.rnorm = rnorm;
    }
    /** 
     * The condition that caused iterations to stop. 
     */
    public Stop stop;
    /**
     * The number of iterations performed.
     */
    public int niter;
    /**
     * The L2 norm ||b|| of the right-hand-side vector b.
     */
    public double bnorm;
    /**
     * The L2 norm ||r|| of the residuals r = b-Ax.
     */
    public double rnorm;
  }

  /**
   * Abstract linear operator A.
   */
  public interface A {

    /**
     * Accumulates the matrix-vector product y = Ax.
     * @param x the input vector x.
     * @param y the output vector y.
     */
    public void apply(Vec x, Vec y);
  }

  /**
   * Constructs a solver with specified parameters.
   * @param tiny threshold for ratio of residuals ||r||/||b||
   * @param maxi maximum number of iterations to perform.
   */
  public CgSolver(double tiny, int maxi) {
    _tiny = tiny;
    _maxi = maxi;
  }

  /**
   * Solves the system of equation Ax = b with CG iterations.
   * @param a the linear operator that represents the matrix A.
   * @param b the right-hand-side vector.
   * @param x the solution vector.
   */
  public Info solve(A a, Vec b, Vec x) {
    return solve(0.0,a,b,x);
  }

  /**
   * Solves the system of equation Ax = b with preconditioned CG iterations.
   * @param a the linear operator that represents the matrix A.
   * @param m the preconditioner that approximates the inverse of A.
   * @param b the right-hand-side vector.
   * @param x the solution vector.
   */
  public Info solve(A a, A m, Vec b, Vec x) {
    return solve(0.0,a,m,b,x);
  }

  /**
   * Solves the system of equation Ax = b with CG iterations.
   * @param anorm estimate for norm ||A|| of linear operator A.
   * @param a the linear operator that represents the matrix A.
   * @param b the right-hand-side vector.
   * @param x the solution vector.
   */
  public Info solve(double anorm, A a, Vec b, Vec x) {
    Vec q = b.clone();
    a.apply(x,q); // q = Ax
    Vec r = b.clone();
    r.add(1.0,q,-1.0); // r = b-Ax
    Vec d = r.clone();
    double bnorm = b.norm2();
    double rnorm = r.norm2();
    double xnorm = x.norm2();
    double rrnorm = rnorm*rnorm;
    logInit(bnorm,rnorm);
    Info info = null;
    int iter;
    for (iter=0; iter<_maxi && rnorm>_tiny*(anorm*xnorm+bnorm); ++iter) {
      logIter(iter,rnorm);
      a.apply(d,q);
      double dq = d.dot(q);
      double alpha = rrnorm/dq;
      x.add(1.0,d,alpha);
      if (anorm>0.0)
        xnorm = x.norm2();
      if (iter%50==49) { // if accumulated rounding error may be large, ...
        a.apply(x,q); // q = Ax
        r.add(0.0,b,1.0); // r = b
        r.add(1.0,q,-1.0); // r = b-Ax
      } else { // otherwise, use shortcut to update residual
        r.add(1.0,q,-alpha); // r -= alpha*q
      }
      double rrnormOld = rrnorm;
      rnorm = r.norm2();
      rrnorm = rnorm*rnorm;
      double beta = rrnorm/rrnormOld;
      d.add(beta,r,1.0);
    }
    logDone(iter, rnorm);
    if (info==null) {
      Stop stop = (iter<_maxi) ? Stop.TINY : Stop.MAXI;
      info = new Info(stop,iter,bnorm,rnorm);
    }
    return info;
  }

  /**
   * Solves the system of equation Ax = b with preconditioned CG iterations.
   * @param anorm estimate for norm ||A|| of linear operator A.
   * @param a the linear operator that represents the matrix A.
   * @param m the preconditioner that approximates the inverse of A.
   * @param b the right-hand-side vector.
   * @param x the solution vector.
   */
  public Info solve(double anorm, A a, A m, Vec b, Vec x) {
    Vec q = b.clone();
    a.apply(x,q); // q = Ax
    Vec r = b.clone();
    r.add(1.0,q,-1.0); // r = r-q = b-Ax
    Vec s = r.clone();
    m.apply(r,s); // s = Mr
    Vec d = s.clone(); // d = s
    double rsnorm = r.dot(s); // r's = r'Mr
    double bnorm = b.norm2();
    double rnorm = r.norm2();
    double xnorm = x.norm2();
    logInit(bnorm,rnorm);
    int iter;
    Info info = null;
    for (iter=0; iter<_maxi && rnorm>_tiny*(anorm*xnorm+bnorm); ++iter) {
      logIter(iter,rnorm);
      a.apply(d, q); // q = Ad
      double dq = d.dot(q); // d'q
      double alpha = rsnorm/dq; // alpha = r'Mr/d'q
      x.add(1.0,d, alpha); // x = x+alpha*d
      xnorm = x.norm2(); // ||x||
      if (iter%50==49) { // if accumulated rounding error may be large, ...
        a.apply(x,q); // q = Ax
        r.add(0.0,b,1.0); // r = b
        r.add(1.0,q,-1.0); // r = b-Ax
      } else { // otherwise, use shortcut to update residual
        r.add(1.0,q,-alpha); // r = r-alpha*q
      }
      rnorm = r.norm2(); // ||r||
      m.apply(r,s); // s = Mr
      double rsnormOld = rsnorm;
      rsnorm = r.dot(s); // r's = r'Mr
      double beta = rsnorm/rsnormOld;
      d.add(beta,s,1.0); // d = s+beta*d
    }
    logDone(iter,rnorm);
    if (info==null) {
      Stop stop = (iter<_maxi) ? Stop.TINY : Stop.MAXI;
      info = new Info(stop,iter,bnorm,rnorm);
    }
    return info;
  }
 
  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _tiny; // converged: norm(r)<tiny*(norm(A)*norm(x)+norm(b))
  private int _maxi; // upper limit on number of iterations

  // Logging.
  private static Logger _log = 
    Logger.getLogger(CgSolver.class.getName());
  private static void logInit(double bnorm, double rnorm) {
    String s = String.format("begin: bnorm=%1.8g rnorm=%1.8g%n",bnorm,rnorm);
    _log.fine(s);
  }
  private static void logIter(int iter, double rnorm) {
    String s = String.format("iter=%d rnorm=%1.8g%n",iter,rnorm);
    _log.finer(s);
  }
  private static void logDone(int iter, double rnorm) {
    String s = String.format("end: iter=%d rnorm=%1.8g%n",iter,rnorm);
    _log.fine(s);
  }
}
