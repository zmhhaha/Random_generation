#ifndef IO_H
#define IO_H
#include <sstream>
#include <iomanip>
#include "geometry.h"
#include "block3D.h"
#include <fstream>


class TeceplotWriter{
public:
    TeceplotWriter(std::string const& fileName_): fileName(fileName_){ 
    }
    ~TeceplotWriter() { }
    template<typename Descriptor>
    void writeData(BlockLattice3D<Descriptor> const& lattice,Box3D const& domain,int offset){
        int Grid_X=domain.getNx();
        int Grid_Y=domain.getNy();
        int Grid_Z=domain.getNz();
        ofstream fout;
        fout.open(fileName+".dat", ios::ate); 
        fout << "TITLE = \"contour\"\nvariables = \"x\", \"y\", \"z\", \"solid\"\nZone I = " << Grid_Z << ", J = " << Grid_Y << ", K = " << Grid_X << " F = POINT" << endl;//按照tecplot的拓扑结构要求输出
        for (int i = domain.x0; i <= domain.x1; i++){
            for (int j = domain.y0; j <= domain.y1; j++){
                for (int k = domain.z0; k <= domain.z1; k++){
                    fout << i << "\t" << j << "\t" << k << "\t" << lattice.get(i, j, k)[offset] << endl;
                }
            }
        }
        fout.close();
    }
    void writeData(ScalarField3D const& lattice,Box3D const& domain){
        int Grid_X=domain.getNx();
        int Grid_Y=domain.getNy();
        int Grid_Z=domain.getNz();
        ofstream fout;
        fout.open(fileName+".dat", ios::ate); 
        fout << "TITLE = \"contour\"\nvariables = \"x\", \"y\", \"z\", \"solid\"\nZone I = " << Grid_Z << ", J = " << Grid_Y << ", K = " << Grid_X << " F = POINT" << endl;//按照tecplot的拓扑结构要求输出
        for (int i = domain.x0; i <= domain.x1; i++){
            for (int j = domain.y0; j <= domain.y1; j++){
                for (int k = domain.z0; k <= domain.z1; k++){
                    fout << i << "\t" << j << "\t" << k << "\t" << lattice.get(i, j, k) << endl;
                }
            }
        }
        fout.close();
    }
private:
    std::string fileName;
};


class PalabosWriter{
public:
    PalabosWriter(std::string const& fileName_): fileName(fileName_){ 
    }
    ~PalabosWriter() { }
    template<typename Descriptor>
    void writeData(BlockLattice3D<Descriptor> const& lattice,Box3D const& domain,int offset){
        ofstream fout;
        fout.open(fileName+".dat", ios::ate);
        for (int i = domain.x0; i <= domain.x1; i++){
            for (int j = domain.y0; j <= domain.y1; j++){
                for (int k = domain.z0; k <= domain.z1; k++){
                    fout << lattice.get(i, j, k)[offset] << endl;
                }
            }
        }
        fout.close();
    }
    void writeData(ScalarField3D const& lattice,Box3D const& domain){
        ofstream fout;
        fout.open(fileName+".dat", ios::ate);
        for (int i = domain.x0; i <= domain.x1; i++){
            for (int j = domain.y0; j <= domain.y1; j++){
                for (int k = domain.z0; k <= domain.z1; k++){
                    fout << lattice.get(i, j, k) << endl;
                }
            }
        }
        fout.close();
    }
private:
    std::string fileName;
};


inline std::string createFileName(std::string name, int number, int width) {
    std::stringstream fNameStream;
    fNameStream << name << std::setfill('0') << std::setw(width) << number;
    return fNameStream.str();
}


#endif //IO_H
