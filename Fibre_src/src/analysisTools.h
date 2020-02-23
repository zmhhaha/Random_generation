#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "block3D.h"
#include "geometry.h"
#include "io.h"
#include "randomnumber.h"
#include <cmath>
#include <limits>
#include <queue>

class Filter{
public:
    Filter(double a,double b,double c):pixel(3,3,3){
        int sum;
        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for (int k = 0; k < 3; k++){
                    sum=0;
                    if(i-1!=0) sum++;
                    if(j-1!=0) sum++;
                    if(k-1!=0) sum++;
                    switch(sum){
                        case 0:
                            pixel.get(i,j,k)=0;
                            break;
                        case 1:
                            pixel.get(i,j,k)=a;
                            break;
                        case 2:
                            pixel.get(i,j,k)=b;
                            break;
                        case 3:
                            pixel.get(i,j,k)=c;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    virtual ~Filter(){}
    double get(int x,int y,int z){
        return pixel.get(x+1,y+1,z+1);
    }
    double const& get(int x,int y,int z) const{
        return pixel.get(x+1,y+1,z+1);
    }
public:
    ScalarField3D pixel;
};

class Node:public Dot3D{
public:
    Node():Dot3D(0,0,0),step(0),parent(NULL){}
    Node(Dot3D const& dot, double step_, Node* parent_=NULL):Dot3D(dot),step(step_),parent(parent_){}
    void setNode(Dot3D const& dot, double step_, Node* parent_){
        x=dot.x;
        y=dot.y;
        z=dot.z;
        step=step_;
        parent=parent_;
    }
    virtual ~Node(){}
public:
    double step;
    Node* parent;
};

enum mode{d6=6,d18=18,d26=26};
class addScale{
public:
    addScale(enum mode m):n(m){}
    std::vector<int> X={0,0,0,0,1,-1,0,0,1,-1,0,0,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,1,1,1};
    std::vector<int> Y={0,0,1,-1,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,-1,-1,1,1,-1,-1,1,1};
    std::vector<int> Z={1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1,0,0,0,0,-1,1,-1,1,-1,1,-1,1};
    int n;
};

class analysisTools
{
public:
    analysisTools(ScalarField3D const& input, double init=0.0, int envelope_=1)
    :envelope(envelope_),
     analysisdomain(input.getBoundingBox().x0,input.getBoundingBox().x1+2*envelope,
                    input.getBoundingBox().y0,input.getBoundingBox().y1+2*envelope,
                    input.getBoundingBox().z0,input.getBoundingBox().z1+2*envelope),
     rawdata(analysisdomain.getNx(),analysisdomain.getNy(),analysisdomain.getNz(),init),
     resultdata(analysisdomain.getNx(),analysisdomain.getNy(),analysisdomain.getNz(),init){
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        for (int i = datadomain.x0; i <= datadomain.x1; i++){
            for (int j = datadomain.y0; j <= datadomain.y1; j++){
                for (int k = datadomain.z0; k <= datadomain.z1; k++){
                    rawdata.get(i,j,k)=input.get(i-envelope, j-envelope, k-envelope);
                }
            }
        }
    }
    ~analysisTools(){}
    void distanceTransform(Filter const& filter){
        resultdata.reset();
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        ScalarField3D first=rawdata;
        resultdata=rawdata;
        int sum,t=0;
        do{
            sum=0;
            t++;
            #pragma omp parallel for reduction(+:sum)
            for (int i = datadomain.x0; i <= datadomain.x1; i++){
                for (int j = datadomain.y0; j <= datadomain.y1; j++){
                    for (int k = datadomain.z0; k <= datadomain.z1; k++){
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
                                        }
                                    }
                                }
                            }
                            if(nsum>0) resultdata.get(i,j,k)=small;
                            else resultdata.get(i,j,k)=0;
                            sum++;
                        }
                    }
                }
            }
            first=resultdata;
        }while(sum!=0);
    }
    double shortestPath(int n, bool printall=false){
        resultdata.reset();
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        RandomNumber rn;
        DotList3D begin;
        ofstream fout;
        fout.open("spdata.dat", ios::ate);
        double sum=0;
        while(begin.getN()<n){
            int tempz=datadomain.z0;
            int tempx=rn.getRandomInt(datadomain.x0,datadomain.x1);
            int tempy=rn.getRandomInt(datadomain.y0,datadomain.y1);
            if(abs(rawdata.get(tempx,tempy,tempz)-0.0)<0.0001){
                begin.addDot(Dot3D(tempx,tempy,tempz));
                double oneans=findShortestPath(Dot3D(tempx,tempy,tempz),printall);
                fout<<oneans<<endl;
                sum+=oneans;
            }
        }
        fout.close();
        return sum/static_cast<double>(n);
    }
    
    void exportAnalysisData(std::string headname, bool palabos=true, bool teceplot=false){
        if(teceplot){
            TeceplotWriter tw("teceplot"+headname);
	        tw.writeData(resultdata,resultdata.getBoundingBox().enlarge(-envelope));
        }
        if(palabos){
            PalabosWriter pw("palabos"+headname);
            pw.writeData(resultdata,resultdata.getBoundingBox().enlarge(-envelope));
        }
    }
private:
    void getPath(Node* parentNode){
        if(!parentNode) return;
        else{
            resultdata.get(parentNode->x,parentNode->y,parentNode->z)=parentNode->step;
            getPath(parentNode->parent);
        }
    }
    double findShortestPath(Dot3D const& begin, bool printall){
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        Node*** parents=new Node**[analysisdomain.getNx()];
        for (int i = 0; i < analysisdomain.getNx(); i++){
            parents[i]=new Node*[analysisdomain.getNy()];
            for (int j = 0; j < analysisdomain.getNy(); j++){
                parents[i][j]=new Node[analysisdomain.getNz()];
            }
        }

        addScale as(mode::d6);

        Filter fl(1,sqrt(2),sqrt(3));
        ScalarField3D inq=rawdata;
        std::queue<Dot3D> q;
        q.push(begin);
        inq.get(begin.x,begin.y,begin.z)=1;
        double length=0;
        while(!q.empty()){
            Dot3D top=q.front();
            q.pop();
            if(top.z!=analysisdomain.z1-1) {
                for(int i=0;i<as.n;i++){
                    int newX=top.x+as.X[i];
                    int newY=top.y+as.Y[i];
                    int newZ=top.z+as.Z[i];
                    if(contained(newX,newY,newZ,datadomain)&&((inq.get(newX,newY,newZ)-0.0)<0.0001)){
                        q.push(Dot3D(newX,newY,newZ));
                        parents[newX][newY][newZ].setNode(Dot3D(newX,newY,newZ),parents[top.x][top.y][top.z].step+fl.get(as.X[i],as.Y[i],as.Z[i]),&parents[top.x][top.y][top.z]);
                        inq.get(newX,newY,newZ)=1;
                    }
                }
            }
            else{
                length=parents[top.x][top.y][top.z].step;
                resultdata.get(top.x,top.y,top.z)=length;
                getPath(parents[top.x][top.y][top.z].parent);
                if(!printall) break;
            }
        }

        for (int i = 0; i < analysisdomain.getNx(); i++){
            for (int j = 0; j < analysisdomain.getNy(); j++){
                delete [] parents[i][j];
            }
        }
        for (int i = 0; i < analysisdomain.getNx(); i++){
            delete [] parents[i];
        }
        delete [] parents;
        return length;
    }
private:
    int envelope;
    Box3D analysisdomain;
    ScalarField3D rawdata;
    ScalarField3D resultdata;
};

#endif //ANALYSISTOOLS_H