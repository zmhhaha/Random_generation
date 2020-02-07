#include "geometry.h"
#include <cstdio>
#include "vec.h"
#include "randomnumber.h"
#include "distribution.h"
#include "cell.h"
#include "descriptor.h"
#include "block3D.h"
#include "fibreset.h"
#include "omp.h"
#include "analysisTools.h"
#include "io.h"

using namespace std;

typedef descriptors::CellDescriptor DESCRIPTOR;

int main(){
    Point3D a(1.25,2.37,3.56),b(3.12,4.34,5.22);
    Dot3D c(1,2,3),d(3,4,5);
    Vec3D v(a,b);
    UnitVec uv;
    a=v.begin;
    b=v.end;
    //a.SetPoint(1,2,3);
    std::cout<<a.x<<std::endl;
    std::cout<<a.y<<std::endl;
    std::cout<<a.z<<std::endl;
    std::cout<<b.x<<std::endl;
    std::cout<<b.y<<std::endl;
    std::cout<<b.z<<std::endl;
    std::cout<<v.norm()<<std::endl;
	std::cout<<uv.norm()<<std::endl;

    RandomNumber randtest;

    std::cout<<randtest.getRandomInt(-500,500)<<std::endl;
    std::cout<<randtest.getRandomZtoO()<<std::endl;

    std::vector<int> vc;
    AnisotropyDistribution ad(20);
    map<double,double> rhomap=ad.getRhoMap();
    map<double,double> dis=ad.getDistribution();
    randtest.getRandomfromDistribution(rhomap);

    ad.exportDistribution();
    ad.exportRhoMap();
	
	BlockLattice3D<DESCRIPTOR> bl(50,50,50);
	BlockLattice3D<DESCRIPTOR> cl(bl);
	bl=cl;

	Cell<DESCRIPTOR> ced;
	Cell<DESCRIPTOR> copy;
	copy.attributeF(ced);

    array<double,2> f=cl.get(2,5,7).getRawParameter();
	double* de=cl.get(2,5,7).getExternal(0);
    printf("%d\t%f\n",f[0],*de);

	Fibre fb(Point3D(1,2,3),Vec3D(Point3D(4,5,6)),2);
	UnitVec direction=fb.getDirection();
	double rbd=direction.norm();
	printf("%f\t%f\t%f\t%f\n",direction.getNx(),direction.getNy(),direction.getNz(),rbd);
	bool infibre=contained(55,45,21,fb);
	printf("%d\n",infibre);

	FibreSet fs;
	fs.insertFibre(fb);

    RandomPoint rp;

    for (int i = 0; i < 100; i++)
    {
        Point3D nrp=rp.NormalRandomPoint(0,100,0,100,0,100,randtest);
        printf("%f\t%f\t%f\n",nrp.x,nrp.y,nrp.z);
        Point3D bp=rp.BetaPoint(1000,rhomap,randtest);
        printf("%f\t%f\t%f\n",bp.x,bp.y,bp.z);
    }
    omp_set_num_threads(36);
    #pragma omp parallel for
    for(int i=0;i<100;i++){
        int id=omp_get_thread_num();
        printf("%d\t",id);
    }
    Filter filter(1,1.414,1.732);
	for(int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
				printf("%d %d %d %f\n",i,j,k,filter.pixel.get(i,j,k));
			}
		}
	}
	

    return 0;
}

/*
int main(int argc, char* argv[])
{
	clock_t starttime, endtime;
	if (argc != 8)
	{
		std::cout << "有一些参数输入错误\n";
		std::cout << "输入结构是：\n";
		std::cout << "1. 纤维半径\n";
		std::cout << "2. 空隙率\n";
		std::cout << "3. 总厚度\n";
		std::cout << "4. 总长度\n";
		std::cout << "5. 总宽度\n";
		std::cout << "6. 精度\n";
		std::cout << "7. 各向异性参数\n";
		exit(EXIT_FAILURE);
	}

	float radius = atof(argv[1]);
	float voidage = atof(argv[2]);
	float thick = atof(argv[3]);
	float length = atof(argv[4]);
	float width = atof(argv[5]);
	float precision = atof(argv[6]);
	float beta = atof(argv[7]);

	starttime = clock();

	int Grid_Z;
	Grid_Z = (int)(thick / precision) + 1;

	int Grid_Y;
	Grid_Y = (int)(length / precision) + 1;

	int Grid_X;
	Grid_X = (int)(width / precision) + 1;

	double radius_grid = radius / precision;

	ofstream fout;

	//打开记录文件
	fout.open("record.dat", ios::ate);
	fout << "纤维半径：" << radius << endl;
	fout << "空隙率：" << voidage << endl;
	fout << "总厚度：" << thick << endl;
	fout << "总长度：" << length << endl;
	fout << "总宽度：" << width << endl;
	fout << "精度：" << precision << endl;
	fout << "Z：" << Grid_Z << endl;
	fout << "Y：" << Grid_Y << endl;
	fout << "X：" << Grid_X << endl;
	fout.close();

	double virtual_long = sqrt(pow(Grid_X, 2) + pow(Grid_Y, 2) + pow(Grid_Z, 2));
	double onefibre = virtual_long*(pow(radius_grid, 2)*3.14159);
	double volume = Grid_X*Grid_Y*Grid_Z*(1 - voidage);
	int t_first = 2*(int)(volume / onefibre);

	BlockLattice3D<DESCRIPTOR> fibreblock(Grid_X,Grid_Y,Grid_Z);
	AnisotropyDistribution dis(beta);
	RandomNumber rng;
	RandomPoint rp;
	FibreSet fs;
	Box3D box=fibreblock.getBoundingBox();
	for (int i = 0; i < t_first; i++){
		/*double px=rng.getRandomfloat(0,Grid_X);
		double py=rng.getRandomfloat(0,Grid_Y);
		double pz=rng.getRandomfloat(0,Grid_Z);
		printf("%f\t%f\t%f\n",px,py,pz);
		Point3D fibrepoint(px,py,pz);
		double vx=rng.getRandomfloat(-1000,1000);
		double vy=rng.getRandomfloat(-1000,1000);
		double vz=tan(rng.getRandomfromDistribution(dis.getRhoMap())-PI/2)*
				  sqrt(pow(vx,2)+pow(vy,2));
		UnitVec fibredirection(Point3D(vx,vy,vz));
		printf("%f\t%f\t%f\n",fibredirection.getNx(),fibredirection.getNy(),fibredirection.getNz());*/

		/*Box3D box(0,Grid_X,0,Grid_Y,0,Grid_Z);
		Point3D fibrepoint=rp.NormalRandomPoint(box.x0,box.x1,box.y0,box.y1,box.z0,box.z1,rng);
		Point3D directionpoint=rp.BetaPoint(1000,dis.getRhoMap(),rng);
		UnitVec fibredirection(directionpoint);
		fs.insertFibre(Fibre(fibrepoint,fibredirection,radius_grid));
		fs.insertFibre(BetaFibre(box,dis.getRhoMap(),1000,radius_grid,rng));
	}
	/*printf("%d\n",fs.getN());
	for(int i = 1; i <= 100; i++){
		printf("%f\t%f\t%f\n",fs[i].getDirection().getNx(),fs[i].getDirection().getNy(),fs[i].getDirection().getNz());
	}
	int sum=0;
	for(int i = 0; i < Grid_X; i++){
		for (int j = 0; j < Grid_Y; j++){
			for(int k = 0; k < Grid_Z; k++){
				Cell<DESCRIPTOR> modify=fibreblock.get(i,j,k);
				for(int f=1;f<=t_first&&modify[0]==0;f++){
					modify[0]=contained(i,j,k,fs[f]);
				}
				fibreblock.get(i,j,k).attributeValues(modify);
				sum+=modify[0];
			}
		}
	}
	double totalblock=Grid_X*Grid_Y*Grid_Z;
	double tempvoidage=1-sum/totalblock;
	printf("%f\n",tempvoidage);
	while(tempvoidage>voidage){
		BetaFibre newfibre(box,dis.getRhoMap(),1000,radius_grid,rng);
		for(int i = 0; i < Grid_X; i++){
			for (int j = 0; j < Grid_Y; j++){
				for(int k = 0; k < Grid_Z; k++){
				
					Cell<DESCRIPTOR> modify=fibreblock.get(i,j,k);
					if(modify[0]==0){
						modify[0]=contained(i,j,k,newfibre);
						fibreblock.get(i,j,k).attributeValues(modify);
						sum+=modify[0];
					}
					tempvoidage=1-(double)sum/totalblock;
					printf("%f\n",tempvoidage);
				}
			}
		}
		fs.insertFibre(newfibre);
	}
	printf("%d\n",fs.getN());
	
	//printf("%f",1-sum/(Grid_X*Grid_Y*Grid_Z));
	fout.open("data.dat", ios::ate); 
	fout << "TITLE = \"contour\"\nvariables = \"x\", \"y\", \"z\", \"solid\"\nZone I = " << Grid_Z << ", J = " << Grid_Y << ", K = " << Grid_X << " F = POINT" << endl;//按照tecplot的拓扑结构要求输出
	for (int i = 0; i < Grid_X; i++)
	{
		for (int j = 0; j < Grid_Y; j++)
		{
			for (int k = 0; k < Grid_Z; k++)
			{
				fout << i << "\t" << j << "\t" << k << "\t" << fibreblock.get(i, j, k)[0] << endl;
			}
		}
	}
	fout.close();
	std::cout << "tecplot3D结构文件写入完毕" << endl;
    return 0;
}*/


typedef descriptors::CellDescriptor DESCRIPTOR;

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		std::cout << "have some mistake in parameter\n";
		std::cout << "the struct of input\n";
		std::cout << "1. fibre set input file\n";
        std::cout << "2. Raw data input file\n";
        std::cout << "3. voidage\n";
        std::cout << "4. beta\n";
		exit(EXIT_FAILURE);
	}
    std::string fNameIn = argv[1];
    std::string fNameInRawData = argv[2];
    float voidage = atof(argv[3]);
    float beta = atof(argv[4]);

    printf("%d\n",1);

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
    ifstream finRaw(fNameInRawData.c_str());
    BlockLattice3D<DESCRIPTOR> rawdata(Grid_X,Grid_Y,Grid_Z);
    for (int i = 0; i < Grid_X; i++){
        for (int j = 0; j < Grid_Y; j++){
            for (int k = 0; k < Grid_Z; k++){
                finRaw >> rawdata.get(i, j, k)[0];
            }
        }
    }

    printf("%d\n",2);

    Filter filter(1,std::sqrt(2),std::sqrt(3));

    printf("%d\n",3);

    ScalarField3D first(Grid_X+2,Grid_Y+2,Grid_Z+2);
    for (int i = 1; i < Grid_X+1; i++){
        for (int j = 1; j < Grid_Y+1; j++){
            for (int k = 1; k < Grid_Z+1; k++){
                first.get(i,j,k)=-rawdata.get(i-1, j-1, k-1)[0];
            }
        }
    }

    printf("%d\n",4);

    ScalarField3D second=first;
    
    int sum,t=0;
    do{
        sum=0;
        t++;
        #pragma omp parallel for reduction(+:sum)
        for (int i = 1; i < Grid_X+1; i++){
            for (int j = 1; j < Grid_Y+1; j++){
                for (int k = 1; k < Grid_Z+1; k++){
                    int nsum=0;
                    double small=std::numeric_limits<double>::max();
                    if(std::abs(first.get(i,j,k)-0.0)<0.0001){
                        for (int ix = -1; ix <= 1; ix++){
                            for (int iy = -1; iy <= 1; iy++){
                                for (int iz = -1; iz <= 1; iz++){
                                    if(std::abs(first.get(i+ix,j+iy,k+iz)-0.0)>0.001) {
                                        nsum++;
                                        double temp=std::max(0.0,first.get(i+ix,j+iy,k+iz))+filter.get(ix,iy,iz);
                                        small=std::min(small,temp);
                                        //big=std::max(big,temp);
                                    }
                                }
                            }
                        }
                        if(nsum>0) second.get(i,j,k)=small;
                        else second.get(i,j,k)=0;
                        //data.setData(first,Box3D(i-1,i+1,j-1,j+1,k-1,k+1));
                        //second.get(i,j,k)=compare(data,filter);
                        //printf("%f\n",second.get(i,j,k));
                        sum++;
                    }
                }
            }
        }
        first=second;
        printf("%d\n",5);
        
    }while(sum!=0);
    TeceplotWriter tw(createFileName("teceplot",1,3));
	tw.writeData(second,second.getBoundingBox());

    PalabosWriter pw(createFileName("palabos",1,3));
	pw.writeData(second,second.getBoundingBox());
    /*data.setData(first,Box3D(5,7,5,7,5,7));

    for (int i = -1; i <=1 ; i++){
        for (int j = -1; j <=1; j++){
            for (int k = -1; k <=1; k++){
                printf("%f\t",data.get(i,j,k));
            }
        }
    }*/
    
    
    printf("%d\n",6);


    /*BlockLattice3D<DESCRIPTOR> surfacedata(Grid_X,Grid_Y,Grid_Z);
    for (int i = 0; i < Grid_X; i++){
        for (int j = 0; j < Grid_Y; j++){
            for (int k = 0; k < Grid_Z; k++){
                finRaw >> surfacedata.get(i, j, k)[1];
            }
        }
    }*/
    //AnisotropyFibreCase<DESCRIPTOR> fibrecase(Grid_X,Grid_Y,Grid_Z,rawdata);
    cout<<"1"<<endl;

	//AnisotropyFibreCase<DESCRIPTOR> fibrecase(Grid_X, Grid_Y, Grid_Z,beta,fs);
	
	/*
	//输出数据tecplot格式

	fout.open("teceplot.dat", ios::ate); 
	fout << "TITLE = \"contour\"\nvariables = \"x\", \"y\", \"z\", \"solid\"\nZone I = " << Grid_Z << ", J = " << Grid_Y << ", K = " << Grid_X << " F = POINT" << endl;//按照tecplot的拓扑结构要求输出
	for (int i = 0; i < Grid_X; i++){
		for (int j = 0; j < Grid_Y; j++){
			for (int k = 0; k < Grid_Z; k++){
				fout << i << "\t" << j << "\t" << k << "\t" << fibrecase.getBlockLattice().get(i, j, k)[1] << endl;
			}
		}
	}
	fout.close();
	//std::cout << "tecplot3D结构文件写入完毕" << endl;*/
	//TeceplotWriter<DESCRIPTOR> tw(createFileName("teceplot",0,0));
	//tw.writeData(fibrecase.getBlockLattice(),fibrecase.getBlockLattice().getBoundingBox(),1);

    //TeceplotWriter<DESCRIPTOR> tw(createFileName("teceplot",0,0));
	//tw.writeData(second,second.getBoundingBox());

	//输出数据palabos格式
	//PalabosWriter<DESCRIPTOR> pw(createFileName("palabos",0,0));
	//pw.writeData(fibrecase.getBlockLattice(),fibrecase.getBlockLattice().getBoundingBox(),1);
	
	/*fout.open("palabos.dat", ios::ate);

	for (int i = 0; i < Grid_X; i++){
		for (int j = 0; j < Grid_Y; j++){
			for (int k = 0; k < Grid_Z; k++){
				fout << fibrecase.getBlockLattice().get(i, j, k)[1] << endl;
			}
		}
	}
	fout.close();*/
	//std::cout << "palabos文件写入完毕" << endl;
    return 0;
}