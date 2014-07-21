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
 * Fault cells in a 3D sampling grid. Each grid sample indexed by (i1,i2,i3)
 * contains either one fault cell or null. The grid facilitates searches for
 * cell nabors in skins and fast iterations along fault traces tangent to
 * fault strike and fault curves tangent to fault dip.
 * <p> 
 * Grid indices need not (and typically do not) begin at zero. Index bounds
 * for a fault cell grid are determined by the minima and maxima of indices of
 * cells used to construct the grid.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.07.06
 */
public class FaultCellGrid {

  /**
   * Constructs a fault grid for specified cells. Grid index bounds are
   * determined by the minimum and maximum indices of the specified cells.
   * @param cells array of cells to be included in the grid.
   */
  public FaultCellGrid(FaultCell[] cells) {
    int i1min = Integer.MAX_VALUE;
    int i2min = Integer.MAX_VALUE;
    int i3min = Integer.MAX_VALUE;
    int i1max = -i1min;
    int i2max = -i2min;
    int i3max = -i3min;
    for (FaultCell cell:cells) {
      if (cell.i1<i1min) i1min = cell.i1;
      if (cell.i2<i2min) i2min = cell.i2;
      if (cell.i3<i3min) i3min = cell.i3;
      if (cell.i1>i1max) i1max = cell.i1;
      if (cell.i2>i2max) i2max = cell.i2;
      if (cell.i3>i3max) i3max = cell.i3;
    }
    _j1 = i1min;
    _j2 = i2min;
    _j3 = i3min;
    _n1 = 1+i1max-i1min;
    _n2 = 1+i2max-i2min;
    _n3 = 1+i3max-i3min;
    _cells = new FaultCell[_n3][_n2][_n1];
    for (FaultCell cell:cells)
      set(cell);
  }

  /**
   * Gets the number of cells in the 1st dimension.
   * @return the number of cells.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of cells in the 2nd dimension.
   * @return the number of cells.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the number of cells in the 3rd dimension.
   * @return the number of cells.
   */
  public int getN3() {
    return _n3;
  }

  /**
   * Gets the lower bound on grid indices in the 1st dimension.
   * @return the lower bound.
   */
  public int getI1Min() {
    return _j1;
  }

  /**
   * Gets the lower bound on grid indices in the 2nd dimension.
   * @return the lower bound.
   */
  public int getI2Min() {
    return _j2;
  }

  /**
   * Gets the lower bound on grid indices in the 3rd dimension.
   * @return the lower bound.
   */
  public int getI3Min() {
    return _j3;
  }

  /**
   * Gets the upper bound on grid indices in the 1st dimension.
   * @return the upper bound.
   */
  public int getI1Max() {
    return _j1+_n1-1;
  }

  /**
   * Gets the upper bound on grid indices in the 2nd dimension.
   * @return the upper bound.
   */
  public int getI2Max() {
    return _j2+_n2-1;
  }

  /**
   * Gets the upper bound on grid indices in the 3rd dimension.
   * @return the upper bound.
   */
  public int getI3Max() {
    return _j3+_n3-1;
  }

  /**
   * Gets the fault cell with specified indices, if any.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @return the fault cell; null, if none or if indices are out of bounds.
   */
  public FaultCell get(int i1, int i2, int i3) {
    i1 -= _j1; 
    i2 -= _j2; 
    i3 -= _j3;
    if (0<=i1 && i1<_n1 && 
        0<=i2 && i2<_n2 && 
        0<=i3 && i3<_n3) {
      return _cells[i3][i2][i1];
    } else {
      return null;
    }
  }

  /**
   * Sets the specified fault cell. Uses the cell's {x1,x2,x3} coordinates to
   * determine the indices of the cell in this grid.
   * @param cell the fault cell.
   */
  public void set(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    _cells[i3][i2][i1] = cell;
  }

  /**
   * Finds a fault cell above the specified cell. Searches for a cell above
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell above,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell above.
   * @return the cell above; null, if none.
   */
  public FaultCell findCellAbove(FaultCell cell) {
    if (cell==null) return null;
    if (cell.ca!=null) return cell.ca;
    return findCellAboveBelow(true,cell);
  }

  /**
   * Finds a fault cell below the specified cell. Searches for a cell below
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell below,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell below.
   * @return the cell below; null, if none.
   */
  public FaultCell findCellBelow(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cb!=null) return cell.cb;
    return findCellAboveBelow(false,cell);
  }

  /**
   * Finds a fault cell left of the specified cell. Searches for a cell left
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * left, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell left.
   * @return the cell left; null, if none.
   */
  public FaultCell findCellLeft(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cl!=null) return cell.cl;
    return findCellLeftRight(true,cell);
  }

  /**
   * Finds a fault cell right of the specified cell. Searches for a cell right
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * right, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell right.
   * @return the cell right; null, if none.
   */
  public FaultCell findCellRight(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cr!=null) return cell.cr;
    return findCellLeftRight(false,cell);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells

  private void init(FaultCell[] cells) {
    if (cells!=null) {
      for (FaultCell cell:cells) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        _cells[i3-_j3][i2-_j2][i1-_j1] = cell;
      }
    }
  }

  private FaultCell findCellAboveBelow(boolean above, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float u1 = cell.u1;
    float u2 = cell.u2;
    float u3 = cell.u3;
    int k1 = 1;
    if (above) {
      k1 = -k1;
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        FaultCell c = get(i1+k1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float du = d1*u1+d2*u2+d3*u3;
          if (du>0.0f) {
            d1 -= du*u1;
            d2 -= du*u2;
            d3 -= du*u3;
            float d = d1*d1+d2*d2+d3*d3; // squared distance to dip line
            if (d<dmin) {
              cmin = c;
              dmin = d;
            }
          }
        }
      }
    }
    return cmin;
  }

  // The search for a cell left or right is not so straightforward as for a
  // cell above or below. The specified cell has eight adjacent samples. We
  // want a cell that is both nearby and located in the strike direction (if
  // right) or opposite direction (if left). We therefore first look for the
  // best cell among the N, E, S, and W adjacent samples, because they are
  // likely to be nearest, and we do not want to skip over them. If and only
  // if we do not find any candidate cells located in the specified direction,
  // we then look among the NE, SE, SW, and NW samples.
  private static final int[] K2LR = { 0, 1, 0,-1, 1, 1,-1,-1};
  private static final int[] K3LR = { 1, 0,-1, 0, 1,-1,-1, 1};
  private FaultCell findCellLeftRight(boolean left, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v1 = cell.v1;
    float v2 = cell.v2;
    float v3 = cell.v3;
    if (left) {
      v1 = -v1;
      v2 = -v2;
      v3 = -v3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int ik=0; ik<8; ++ik) {
      if (ik==4 && cmin!=null)
        break;
      int k2 = K2LR[ik];
      int k3 = K3LR[ik];
      FaultCell c = get(i1,i2+k2,i3+k3);
      if (c!=null) {
        float d1 = c.x1-x1;
        float d2 = c.x2-x2;
        float d3 = c.x3-x3;
        float dv = d1*v1+d2*v2+d3*v3;
        if (dv>0.0f) {
          d1 -= dv*v1;
          d2 -= dv*v2;
          d3 -= dv*v3;
          float d = d1*d1+d2*d2+d3*d3; // squared distance to strike line
          if (d<dmin) {
            cmin = c;
            dmin = d;
          }
        }
      }
    }
    return cmin;
  }
}
