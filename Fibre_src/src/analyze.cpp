#include <cstdio>
#include <fstream>
#include <ctime>
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
	if (argc != 7)
	{
		std::cout << "have some mistake in parameter\n";
		std::cout << "the struct of input\n";
		std::cout << "1. fibre set input file\n";
        std::cout << "2. Raw data input file\n";
        std::cout << "3. voidage\n";
        std::cout << "4. beta\n";
        std::cout << "5. analytical model\n";
        std::cout << "6. the other parameter\n";
		exit(EXIT_FAILURE);
	}

    std::string fNameIn = argv[1];
    std::string fNameInRawData = argv[2];
    float voidage = atof(argv[3]);
    float beta = atof(argv[4]);
    int mod = atoi(argv[5]);
    int typeparameter = atoi(argv[6]);
    std::string inputfile = argv[6];

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
    fin.close();

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
    finRaw.close();

    Box3D analyzedDomain(0,Grid_X-1,0,Grid_Y-1,0,Grid_Z-1);
    analysisTools astools(*extractSubDomain(rawdata,analyzedDomain));
    printf("Raw Data input.\n");

    Filter filter(1,std::sqrt(2),std::sqrt(3));
    

    printf("Begin analyze.\n");
    if(mod==1){
        addScale addscale(d26);
        astools.distanceTransform(filter,addscale,2);
        outname="df";
    }
    if(mod==2){
        int n=1000;
        double nsp=astools.shortestPath(n);
        printf("Statistical sample size:%d\n",n);
        printf("The shortest Path:%f\n",nsp);
        outname="sp";
    }
    if(mod==3){
        int n=1000;
        int length=typeparameter;
        double nsp=astools.porousDistribution(n,length);
        printf("Statistical sample size:%d\n",n);
        printf("The average voidage:%f\n",nsp);
        return 0;
    }
    if(mod==4){
        int part=20;
        addScale addscale(d6);
        int nc=astools.porousCenter(part,filter,addscale,0);
        printf("Block number:%d\n",part);
        printf("The center number:%d\n",nc);
        return 0;
    }
    if(mod==5){
        time_t time_start=time(NULL);
        addScale addscale(d26);
        double nc=astools.smallWeightPath(filter,addscale,typeparameter);
        printf("The average path length:%f\n",nc);
        outname="swp";
        time_t time_end=time(NULL);
        cout<<"time use:"<<(double)difftime(time_end,time_start)<<"s"<<endl;
    }
    if(mod==6){
        bool printbool=true;
        ifstream fin(inputfile.c_str());
        string filename=inputfile.substr(0,inputfile.size()-4);
        int thetempvalue;
        vector<int> pointvalue;
        while(fin>>thetempvalue){
            pointvalue.push_back(thetempvalue);
        }
        fin.close();
        DotList3D pointlist;
        for(int i=0;i<pointvalue.size();i+=3){
            pointlist.addDot(Dot3D(pointvalue[i],pointvalue[i+1],pointvalue[i+2]));
        }
        addScale addscale(d26);
        astools.printdistancevalue(filter, addscale, pointlist,filename,printbool);
        if(!printbool) return 0;
        outname="pcball";
        outname+=filename;
    }

    astools.exportAnalysisData(outname);
    return 0;
}