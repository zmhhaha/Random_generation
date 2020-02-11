#ifndef ANISOTROPYFIBRECASE_H
#define ANISOTROPYFIBRECASE_H

#include "geometry.h"
#include <cstdio>
#include <fstream>
#include "vec.h"
#include "randomnumber.h"
#include "distribution.h"
#include "cell.h"
#include "descriptor.h"
#include "block3D.h"
#include "fibreset.h"
#include "omp.h"

template <typename Descriptor>
class AnisotropyFibreCase
{
public:
    AnisotropyFibreCase(int Grid_X, int Grid_Y, int Grid_Z,double voidage, double Grid_Radius, double beta_)
                        :fibreblock(Grid_X,Grid_Y,Grid_Z),dis(beta_){
        setFibreSet(voidage,Grid_Radius);
        printFibreBlock();
        int sum=checkSolid();
        double tempvoidage=1-(double)sum/fibreblock.getBoundingBox().nCells();
        printf("%f\n",tempvoidage);
        amendFibreSet(voidage,Grid_Radius);
        printf("%d\n",fs.getN());
        checkSurface();
        exportFibreSetParameter();
        exportRawData();
    }
    AnisotropyFibreCase(int Grid_X, int Grid_Y, int Grid_Z,double beta_,FibreSet fs_)
                        :fibreblock(Grid_X,Grid_Y,Grid_Z),dis(beta_),fs(fs_){
        printFibreBlock();
        checkSurface();
        exportRawData();
    }
    AnisotropyFibreCase(int Grid_X, int Grid_Y, int Grid_Z,BlockLattice3D<Descriptor> lattice_):fibreblock(lattice_){
        checkSurface();
    }
    ~AnisotropyFibreCase() { }
    void printFibreBlock() {
        int Grid_X=fibreblock.getBoundingBox().getNx();
        int Grid_Y=fibreblock.getBoundingBox().getNy();
        int Grid_Z=fibreblock.getBoundingBox().getNz();
        #pragma omp parallel for
        for(int i = 0; i < Grid_X; i++){
            for (int j = 0; j < Grid_Y; j++){
                for(int k = 0; k < Grid_Z; k++){
                    Cell<Descriptor> modify=fibreblock.get(i,j,k);
                    for(int f=1;f<=fs.getN()&&modify[0]==0;f++){
                        modify[0]=contained(i,j,k,fs[f]);
                    }
                    fibreblock.get(i,j,k).attributeValues(modify);
                }
            }
	    }
    }
    void amendFibreSet(double voidage,double Grid_Radius){
        int Grid_X=fibreblock.getBoundingBox().getNx();
        int Grid_Y=fibreblock.getBoundingBox().getNy();
        int Grid_Z=fibreblock.getBoundingBox().getNz();
        int sum=checkSolid();
        int totalN=fibreblock.getBoundingBox().nCells();
        double tempvoidage=1-(double)sum/(double)totalN;

        while(tempvoidage>voidage){
            BetaFibre newfibre(fibreblock.getBoundingBox(),dis.getRhoMap(),1000,Grid_Radius,rng);
            int tempsum=0;
            #pragma omp parallel for reduction(+:tempsum)
            for(int i = 0; i < Grid_X; i++){
                for (int j = 0; j < Grid_Y; j++){
                    for(int k = 0; k < Grid_Z; k++){
                        Cell<Descriptor> modify=fibreblock.get(i,j,k);
                        if(modify[0]==0){
                            modify[0]=contained(i,j,k,newfibre);
                            fibreblock.get(i,j,k).attributeValues(modify);
                            tempsum+=modify[0];
                        }
                    }
                }
            }
            sum+=tempsum;
            tempvoidage=1-(double)sum/(double)totalN;
            printf("%f\n",tempvoidage);
            fs.insertFibre(newfibre);
        }
    }
    void setFibreSet(double voidage,double Grid_Radius) {
        int Grid_X=fibreblock.getBoundingBox().getNx();
        int Grid_Y=fibreblock.getBoundingBox().getNy();
        int Grid_Z=fibreblock.getBoundingBox().getNz();
        double virtual_long = sqrt(pow(Grid_X, 2) + pow(Grid_Y, 2) + pow(Grid_Z, 2));
	    double onefibre = virtual_long*(pow(Grid_Radius, 2)*3.14159);
	    double volume = Grid_X*Grid_Y*Grid_Z*(1 - voidage);
	    int t_first = 2*(int)(volume / onefibre);
        Box3D box=fibreblock.getBoundingBox();
        for (int i = 0; i < t_first; i++){
            fs.insertFibre(BetaFibre(box,dis.getRhoMap(),1000,Grid_Radius,rng));
        }
    }
    BlockLattice3D<Descriptor> const& getBlockLattice() const{
        return fibreblock;
    }
    AnisotropyDistribution const& getDistribution() const{
        return dis;
    }
    FibreSet const& getFibreSet() const{
        return fs;
    }
    void checkSurface(){
        int Grid_X=fibreblock.getBoundingBox().getNx();
        int Grid_Y=fibreblock.getBoundingBox().getNy();
        int Grid_Z=fibreblock.getBoundingBox().getNz();
        for (int i = 0; i < Grid_X; i++){
            for (int j = 0; j < Grid_Y; j++){
                for (int k = 0; k < Grid_Z; k++){
                    if (fibreblock.get(i,j,k)[0]==1){
                        int sum=0;
                        for (int ix = -1; ix <= 1; ix++){
                            for (int iy = -1; iy <= 1; iy++){
                                for (int iz = -1; iz <= 1; iz++){
                                    if(contained(i+ix,j+iy,k+iz,fibreblock.getBoundingBox())) sum+=fibreblock.get(i+ix,j+iy,k+iz)[0];
                                }
                            }
                        }
                        if(sum==27) fibreblock.get(i,j,k)[1]=2;
                        else fibreblock.get(i,j,k)[1]=1;
                    }
                }   
            }
        }
    }
    int checkSolid(){
        int Grid_X=fibreblock.getBoundingBox().getNx();
        int Grid_Y=fibreblock.getBoundingBox().getNy();
        int Grid_Z=fibreblock.getBoundingBox().getNz();
        int sum=0;
        #pragma omp parallel for reduction(+:sum)
        for(int i = 0; i < Grid_X; i++){
            for (int j = 0; j < Grid_Y; j++){
                for(int k = 0; k < Grid_Z; k++){
                    sum+=fibreblock.get(i,j,k)[0];
                }
            }
	    }
        return sum;
    }
    void exportFibreSetParameter(){
        ofstream fout;
        fout.open("FibreParameter.dat", ios::ate);
        fout<<fibreblock.getNx()<<"\t"<<fibreblock.getNy()<<"\t"<<fibreblock.getNz()<<endl;
        fout<<fs.getN()<<endl;
        for ( int f = 1; f <= fs.getN(); f++){
            fout<<fs[f].getBasePoint().x<<"\t"<<fs[f].getBasePoint().y<<"\t"<<fs[f].getBasePoint().z<<"\t"
                <<fs[f].getDirection().getNx()<<"\t"<<fs[f].getDirection().getNy()<<"\t"<<fs[f].getDirection().getNz()<<"\t"
                <<fs[f].getRadius()<<endl;
        }
        
        fout.close();
    }
    void exportRawData(){
        int Grid_X=fibreblock.getNx();
        int Grid_Y=fibreblock.getNy();
        int Grid_Z=fibreblock.getNz();
        ofstream fout;
        fout.open("RawData.dat", ios::ate); 
        for (int i = 0; i < Grid_X; i++){
            for (int j = 0; j < Grid_Y; j++){
                for (int k = 0; k < Grid_Z; k++){
                    fout << fibreblock.get(i, j, k)[0] << "\t" << endl;
                }
            }
        }
    }
private:
    BlockLattice3D<Descriptor> fibreblock;
    AnisotropyDistribution dis;
    RandomNumber rng;
    FibreSet fs;
};


#endif //ANISOTROPYFIBRECASE_H