#ifndef BLOCK_3D_H
#define BLOCK_3D_H

#include "geometry.h"
#include "cell.h"


class Block3D {
public:
    Block3D(int nx_, int ny_, int nz_): nx(nx_), ny(ny_), nz(nz_){ }
    virtual ~Block3D() { }
    virtual Box3D getBoundingBox() const{return Box3D(0, getNx()-1, 0, getNy()-1, 0, getNz()-1);}
    int getNx() const { return nx; }
    /// Get number of cells in y-direction.
    int getNy() const { return ny; }
    /// Get number of cells in z-direction.
    int getNz() const { return nz; }
    void swap(Block3D& rhs){
        std::swap(nx, rhs.nx);
        std::swap(ny, rhs.ny);
        std::swap(nz, rhs.nz);
    }
private:
    int nx, ny, nz;/// Dimensions of the lattice.
};

template <typename Descriptor>
class BlockLatticeBase3D {
public:
    BlockLatticeBase3D(){}
    virtual ~BlockLatticeBase3D(){}
public:
    virtual Cell<Descriptor>& get(int iX, int iY, int iZ) =0;
    virtual Cell<Descriptor> const& get(int iX, int iY, int iZ) const =0;
};

class ScalarFieldBase3D {
public:
    virtual ~ScalarFieldBase3D() { }
public:
    virtual void reset() =0;
    virtual double& get(int iX, int iY, int iZ) =0;
    virtual double const& get(int iX, int iY, int iZ) const =0;
};

template <typename Descriptor>
class BlockLattice3D:public BlockLatticeBase3D<Descriptor>,public Block3D
{
public:
    /// Construction of an nx_ by ny_ by nz_ lattice
    BlockLattice3D(int nx_, int ny_, int nz_):Block3D(nx_,ny_,nz_){
        // Allocate memory, and initialize dynamics.
        allocateAndInitialize();
    };
    /// Destruction of the lattice
    virtual ~BlockLattice3D(){
        releaseMemory();
    }
    /// Copy construction
    BlockLattice3D(BlockLattice3D<Descriptor> const& rhs):Block3D(rhs){
        int nx = this->getNx();
        int ny = this->getNy();
        int nz = this->getNz();
        allocateAndInitialize();
        for (int iX=0; iX<nx; ++iX) {
            for (int iY=0; iY<ny; ++iY) {
                for (int iZ=0; iZ<nz; ++iZ) {
                    Cell<Descriptor>& cell = grid[iX][iY][iZ];
                    // Assign cell from rhs
                    cell = rhs.grid[iX][iY][iZ];
                }
            }
        }
    }
    /// Copy assignment
    BlockLattice3D& operator=(BlockLattice3D<Descriptor> const& rhs){
        BlockLattice3D<Descriptor> tmp(rhs);
        swap(tmp);
        return *this;
    }
public:
    /// Read/write access to lattice cells
    virtual Cell<Descriptor>& get(int iX, int iY, int iZ) {
        return grid[iX][iY][iZ];
    }
    /// Read only access to lattice cells
    virtual Cell<Descriptor> const& get(int iX, int iY, int iZ) const {
        return grid[iX][iY][iZ];
    }
    void swap(BlockLattice3D& rhs) {
        Block3D::swap(rhs);
        std::swap(rawData, rhs.rawData);
        std::swap(grid, rhs.grid);
    }
private:
    /// Helper method for memory allocation
    void allocateAndInitialize(){
        int nx = this->getNx();
        int ny = this->getNy();
        int nz = this->getNz();
        rawData = new Cell<Descriptor> [nx*ny*nz];
        grid    = new Cell<Descriptor>** [nx];
        for (int iX=0; iX<nx; ++iX) {
            grid[iX] = new Cell<Descriptor>* [ny];
            for (int iY=0; iY<ny; ++iY) {
                grid[iX][iY] = rawData + nz*(iY+ny*iX);
            }
        }
    }
    /// Helper method for memory de-allocation
    void releaseMemory(){
        int nx = this->getNx();
        int ny = this->getNy();
        int nz = this->getNz();
        delete [] rawData;
        for (int iX=0; iX<nx; ++iX) {
            delete [] grid[iX];
        }
        delete [] grid;
    }
private:
    Cell<Descriptor>     *rawData;
    Cell<Descriptor>   ***grid;
};

class ScalarField3D : public ScalarFieldBase3D, public Block3D {
public:
    ScalarField3D(int nx_, int ny_, int nz_,double iniVal=0):Block3D(nx_, ny_, nz_){
        allocateMemory();
        for (unsigned int iData=0; iData<getSize(); ++iData) {
            (*this)[iData] = iniVal;
        }
    }
    ~ScalarField3D(){
        releaseMemory();
    }
    ScalarField3D(ScalarField3D const& rhs):Block3D(rhs){
        allocateMemory();
        for (unsigned int iData=0; iData<getSize(); ++iData) {
            (*this)[iData] = rhs[iData];
        }
    }
    ScalarField3D& operator=(ScalarField3D const& rhs){
        ScalarField3D tmp(rhs);
        swap(tmp);
        return *this;
    }
    void swap(ScalarField3D& rhs){
        Block3D::swap(rhs);
        std::swap(rawData, rhs.rawData);
        std::swap(field, rhs.field);
    }
public:
    virtual void reset(){
        for (unsigned int index=0; index<this->getNx()*this->getNy()*this->getNz(); ++index) {
            (*this)[index] = 0;
        }
    }
    virtual void reset(double value){
        for (unsigned int index=0; index<this->getNx()*this->getNy()*this->getNz(); ++index) {
            (*this)[index] = value;
        }
    }
    virtual unsigned int getSize() const { return (unsigned int)this->getNx()*(unsigned int)this->getNy()*(unsigned int)this->getNz(); }
    virtual double& get(int iX, int iY, int iZ) {
        return field[iX][iY][iZ];
    }
    virtual double const& get(int iX, int iY, int iZ) const {
        return field[iX][iY][iZ];
    }
    double& operator[] (int ind) {
        return rawData[ind];
    }
    double const& operator[] (int ind) const {
        return rawData[ind];
    }
private:
    void allocateMemory(){
        rawData = new double [(unsigned int)this->getNx()*(unsigned int)this->getNy()*(unsigned int)this->getNz()];
        field   = new double** [(unsigned int)this->getNx()];
        for (int iX=0; iX<this->getNx(); ++iX) {
            field[iX] = new double* [(unsigned int)this->getNy()];
            for (int iY=0; iY<this->getNy(); ++iY) {
                field[iX][iY] = rawData + (unsigned int)this->getNz()*((unsigned int)iY+(unsigned int)this->getNy()*(unsigned int)iX);
            }
        }
    }
    void releaseMemory(){
        for (int iX=0; iX<this->getNx(); ++iX) {
            delete [] field[iX];
        }
        delete [] field;
        delete [] rawData; rawData = 0;
    }
private:
    double   *rawData;
    double   ***field;
};

#endif