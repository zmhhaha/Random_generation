#include <cstdio>
#include <fstream>
#include "geometry.h"
#include "vec.h"
#include "randomnumber.h"
#include "distribution.h"
#include "cell.h"
#include "descriptor.h"
#include "block3D.h"
#include "fibreset.h"
#include "AnisotropyFibreCase.h"
#include "io.h"

#include "analysisTools.h"

using namespace std;

typedef descriptors::CellDescriptor DESCRIPTOR;

int main(int argc, char* argv[])
{
	if (argc != 6)
	{
		std::cout << "have some mistake in parameter\n";
		std::cout << "the struct of input\n";
		std::cout << "1. fibre set input file\n";
        std::cout << "2. Raw data input file\n";
        std::cout << "3. voidage\n";
        std::cout << "4. beta\n";
        std::cout << "5. analytical model\n";
		exit(EXIT_FAILURE);
	}

    std::string fNameIn = argv[1];
    std::string fNameInRawData = argv[2];
    float voidage = atof(argv[3]);
    float beta = atof(argv[4]);
    int mod = atoi(argv[5]);

    std::string outname;

    printf("Parameter input.\n");

    ifstream fin(fNameIn.c_str());
    int Grid_X, Grid_Y, Grid_Z;
    int fibreNumber;
    FibreSet fs;
    fin>>Grid_X>>Grid_Y>>Grid_Z;
    fin>>fibreNumber;
    for (int i = 0; i < fibreNumber; i++){
        double px,py,pz,vx,vy,vz,radius;
        fin>>px>>py>>pz>>vx>>vy>>vz>>radius;
        fs.insertFibre(Fibre(Point3D(px,py,pz),UnitVec(Point3D(vx,vy,vz)),radius));
    }

    printf("Fibre set input.\n");

    ifstream finRaw(fNameInRawData.c_str());
    ScalarField3D rawdata(Grid_X,Grid_Y,Grid_Z);
    for (int i = 0; i < Grid_X; i++){
        for (int j = 0; j < Grid_Y; j++){
            for (int k = 0; k < Grid_Z; k++){
                finRaw >> rawdata.get(i, j, k);
            }
        }
    }

    analysisTools astools(rawdata);
    printf("Raw Data input.\n");

    Filter filter(1,std::sqrt(2),std::sqrt(3));

    printf("Begin analyze.\n");
    if(mod==1){
        astools.distanceTransform(filter);
        outname="df";
    }
    if(mod==2){
        int n=1000;
        double nsp=astools.shortestPath(1000);
        printf("Statistical sample size:%d\n",n);
        printf("The shortest Path:%f\n",nsp);
        outname="sp";
    }

    astools.exportAnalysisData(outname);
    return 0;
}