#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include<vector>
#include<cmath>
#include"template.h"
#include"vec.h"

using std::pow;

class Box3D {
public:
    Box3D() : x0(), x1(), y0(), y1(), z0(), z1() { }
    Box3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
        : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
    { }

    /// Add a border of nCells cells to the box
    Box3D enlarge(int nCells) const {
        return Box3D(x0-nCells, x1+nCells, y0-nCells, y1+nCells, z0-nCells, z1+nCells);
    }

    /// Add a border of nCells cells to the box in one direction (x=0, y=1, z=2)
    Box3D enlarge(int nCells, int dir) const {
        switch ( dir ) {
            case 0:
              return Box3D(x0-nCells, x1+nCells, y0, y1, z0, z1);
            case 1:
              return Box3D(x0, x1, y0-nCells, y1+nCells, z0, z1);
            case 2:
              return Box3D(x0, x1, y0, y1, z0-nCells, z1+nCells);
            default:
              break;
        }
        return Box3D(-1,-1,-1,-1,-1,-1);
    }

    /// Number of cells in x-direction
    int getNx()  const { return (x1-x0+1); }
    /// Number of cells in y-direction
    int getNy()  const { return (y1-y0+1); }
    /// Number of cells in z-direction
    int getNz()  const { return (z1-z0+1); }
    /// Total number of cells in the box
    int nCells() const { return getNx()*getNy()*getNz(); }
    /// Return the maximum of getNx(), getNy(), and getNz()
    int getMaxWidth() const { return std::max(std::max(getNx(), getNy()), getNz()); }
    /// Copy the data into a 6-element array.

    bool operator==(Box3D const& rhs) const {
        return x0 == rhs.x0 && y0 == rhs.y0 && z0 == rhs.z0 &&
               x1 == rhs.x1 && y1 == rhs.y1 && z1 == rhs.z1;
    }

    int x0, x1, y0, y1, z0, z1;
};

/// Coordinates of a 3D point
typedef PointTemplate3D<int> Dot3D;

/// List of 3D points, used to describe a subdomain
class DotList3D {
public:
    DotList3D() { }
    DotList3D(std::vector<Dot3D> const& dots_) : dots(dots_) { }
    /// Add one more point to the list
    void addDot(Dot3D dot) {
        dots.push_back(dot);
    }
    /// Erase the end dot
    void eraseDot() {
        dots.pop_back();
    }
    /// Add more points to the list
    void addDots(std::vector<Dot3D> const& dots_) {
        dots.insert(dots.end(), dots_.begin(), dots_.end());
    }
    /// Get const reference to one of the dots
    Dot3D const& getDot(int whichDot) const {
        return dots[whichDot];
    }
    /// Get non-const reference to one of the dots
    Dot3D& getDot(int whichDot) {
        return dots[whichDot];
    }

    /// Get total number of points
    int getN() const {
        return dots.size();
    }

    std::vector<Dot3D> dots;
};

/// Decide if lattice point is contained in 3D box, boundaries inclusive
inline bool contained(int x, int y, int z, Box3D const& box) {
    return x>=box.x0 && x<=box.x1 &&
           y>=box.y0 && y<=box.y1 &&
           z>=box.z0 && z<=box.z1;
}



#endif