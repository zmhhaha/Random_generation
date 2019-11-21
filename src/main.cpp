#include "geometry.h"
#include <cstdio>
#include "vec.h"
#include "randomnumber.h"
#include "distribution.h"
#include "cell.h"
#include "descriptor.h"
#include "block3D.h"
#include "fibreset.h"
#include "AnisotropyFibreCase.h"

using namespace std;

typedef descriptors::CellDescriptor DESCRIPTOR;

int main(int argc, char* argv[])
{
	clock_t starttime, endtime;
	if (argc != 8)
	{
		std::cout << "have some mistake in parameter\n";
		std::cout << "the struct of input\n";
		std::cout << "1. fibre radius\n";
		std::cout << "2. voidage\n";
		std::cout << "3. total thickness\n";
		std::cout << "4. total length\n";
		std::cout << "5. total width\n";
		std::cout << "6. precision\n";
		std::cout << "7. anisotropy parameter\n";
		std::cout << "the standard input: ./main 1 0.6 10 20 30 0.2 20\n";
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

	AnisotropyFibreCase<DESCRIPTOR> fibrecase(Grid_X, Grid_Y, Grid_Z,voidage, radius_grid, beta);
	
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
	std::cout << "tecplot3D结构文件写入完毕" << endl;

	//输出数据palabos格式
	fout.open("palabos.dat", ios::app); //ios::app 所有输出附加在文件末尾

	for (int i = 0; i < Grid_X; i++){
		for (int j = 0; j < Grid_Y; j++){
			for (int k = 0; k < Grid_Z; k++){
				fout << fibrecase.getBlockLattice().get(i, j, k)[1] << endl;
			}
		}
	}
	fout.close();
	std::cout << "palabos文件写入完毕" << endl;
	/*
	for (int i = 0; i < Grid_X; i++){
		for (int j = 0; j < Grid_Y; j++){
			for (int k = 0; k < Grid_Z; k++){
				cout << i << "\t" << j << "\t" << k << "\t" << fibrecase.getBlockLattice().get(i, j, k)[1] << endl;
			}
		}
	}
	*/
    return 0;
}