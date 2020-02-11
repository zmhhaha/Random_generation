#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#include<fstream>
#include<iostream>
#include<map>
using std::map;
using std::pair;
using std::string;
using std::ofstream;
using std::ios;
using std::endl;
using std::cout;

class Distribution{
protected:
    map<double,double> rhomap;
    map<double,double> dis;
public:
    virtual ~Distribution(){}
    virtual void init()=0;
    map<double,double> const& getRhoMap() const{
        return rhomap;
    }
    map<double,double> const& getDistribution() const{
        return dis;
    }
    virtual void exportRhoMap()=0;
    virtual void exportDistribution()=0;
};

class AnisotropyDistribution:public Distribution
{
private:
    double beta;
public:
    AnisotropyDistribution(double beta_):beta(beta_){
        this->init();
    }
    AnisotropyDistribution():beta(1){
        this->init();
    }
    virtual ~AnisotropyDistribution(){}
    void init();
    void exportRhoMap();
    void exportDistribution();
};

void AnisotropyDistribution::init(){
    int n = 10000;
    double dtheta = PI / n;
    double theta;
    double costheta;
	double sintheta;
    double rho;
    double sumrho=0;
    for (int i = 0; i <= n; i++)
    {
        theta=dtheta*i;
        costheta = cos(theta);
        sintheta = sin(theta);
        rho= 1 * beta * sintheta / (2 * pow((1 + (pow(beta, 2) - 1)*pow(costheta, 2)), 1.5));
        sumrho+=rho*dtheta;
        dis.insert(pair<double,double>(theta,rho));
        rhomap.insert(pair<double,double>(sumrho,theta));
    }
}

void AnisotropyDistribution::exportRhoMap(){
    int n=rhomap.size();
    ofstream fout;
    fout.open("rhomap.dat", ios::ate); 
    fout << "TITLE = \"contour\"\nvariables = \"theta\", \"rhosum\"\nZone I = " << n << " F = POINT" <<endl;//按照tecplot的拓扑结构要求输出
    for (auto c:rhomap)
    {
        fout << c.second << "\t" << c.first << endl;
    }
    fout.close();
}
    
void AnisotropyDistribution::exportDistribution(){
    int n=dis.size();
    ofstream fout;
    fout.open("distribution.dat", ios::ate); 
    fout << "TITLE = \"contour\"\nvariables = \"theta\", \"rho\"\nZone I = " << n << " F = POINT" <<endl;//按照tecplot的拓扑结构要求输出
    for (auto c:dis)
    {
        fout << c.first << "\t" << c.second << endl;
    }
    fout.close();
}

/*
//工厂模式探索
template<typename TheDistribution>
class DistributionCreator{
public:
    DistributionCreator(){}
    virtual ~DistributionCreator(){}
    virtual Distribution* CreateProduct(){
        return new TheDistribution;
    }
};
*/

#endif //DISTRIBUTION_H