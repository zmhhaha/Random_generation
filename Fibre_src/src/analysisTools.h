#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "block3D.h"
#include "geometry.h"
#include "io.h"
#include "randomnumber.h"
#include <cmath>
#include<cfloat>
#include <limits>
#include <queue>
#include <vector>
#include <sstream>
using std::vector;

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

    void distanceTransform(Filter const& filter, addScale const& addscale, int rawdataType){
        resultdata.reset();
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        ScalarField3D first(analysisdomain.getNx(),analysisdomain.getNy(),analysisdomain.getNz(),-1);
        
        // the change of rawdata input;
        // value 0 is the datadomain without wall
        // value 1 is the datadomain with all wall
        // value other is the datadomain with x-directcion and y-direction wall 
        Box3D rawdatainput;
        if(rawdataType==0) rawdatainput=analysisdomain;
        else if(rawdataType==1) rawdatainput=datadomain;
        else rawdatainput=Box3D(datadomain.x0,datadomain.x1,datadomain.y0,datadomain.y1,analysisdomain.z0,analysisdomain.z1);
        //printf("%d,%d,%d,%d,%d,%d\n",rawdatainput.x0,rawdatainput.x1,rawdatainput.y0,rawdatainput.y1,rawdatainput.z0,rawdatainput.z1);
        
        #pragma omp parallel for
        for (int i = rawdatainput.x0; i <= rawdatainput.x1; i++){
            for (int j = rawdatainput.y0; j <= rawdatainput.y1; j++){
                for (int k = rawdatainput.z0; k <= rawdatainput.z1; k++){
                    first.get(i,j,k)=-rawdata.get(i,j,k);
                }
            }
        }
        resultdata=first;
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
                            /*for (int ix = -1; ix <= 1; ix++){
                                for (int iy = -1; iy <= 1; iy++){
                                    for (int iz = -1; iz <= 1; iz++){
                                        if(std::abs(first.get(i+ix,j+iy,k+iz)-0.0)>0.0001) {
                                            nsum++;
                                            double temp=std::max(0.0,first.get(i+ix,j+iy,k+iz))+filter.get(ix,iy,iz);
                                            small=std::min(small,temp);
                                        }
                                    }
                                }
                            }*/
                            for(int as=0;as<addscale.n;as++){
                                int ix=addscale.X[as];
                                int iy=addscale.Y[as];
                                int iz=addscale.Z[as];
                                if(std::abs(first.get(i+ix,j+iy,k+iz)-0.0)>0.0001) {

                                    nsum++;
                                    double temp=std::max(0.0,first.get(i+ix,j+iy,k+iz))+filter.get(ix,iy,iz);
                                    small=std::min(small,temp);
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

        #pragma omp parallel for
        for (int i = datadomain.x0; i <= datadomain.x1; i++){
            for (int j = datadomain.y0; j <= datadomain.y1; j++){
                for (int k = datadomain.z0; k <= datadomain.z1; k++){
                    if(resultdata.get(i,j,k)<0) resultdata.get(i,j,k)=0;
                }
            }
        }
        return ;
    }

    int porousCenter(int part, Filter const& filter, addScale const& addScale,int porousCenterType){
        distanceTransform(filter,addScale,porousCenterType);
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        
        //the compare domain size
        int compareLength=1;

        int Nx=datadomain.getNx();
        int Ny=datadomain.getNy();
        int Nz=datadomain.getNz();
        int interpolation=Nx/part;
        ofstream fout;
        fout.open("pcdata.dat", ios::ate);
        int allcenter=0;
        for(int p=0;p<part-1;p++){
            int npx0 = datadomain.x0+p*interpolation;
            int npx1 = datadomain.x0+(p+1)*interpolation-1;
            std::vector<Dot3D> buffer(interpolation*Ny*Nz,Dot3D(0,0,0));
            std::vector<bool> bufferIndex(interpolation*Ny*Nz,false);
            int numpoint=0;
            #pragma omp parallel for reduction(+:numpoint)
            for(int i=0;i<interpolation;i++){
                for(int j=0;j<Ny;j++){
                    for(int k=0;k<Nz;k++){
                        int index=i*Ny*Nz+j*Nz+k;
                        if(maxpoint(Dot3D(npx0+i,datadomain.y0+j,datadomain.z0+k),datadomain,compareLength)){
                            buffer[index]=Dot3D(npx0+i,datadomain.y0+j,datadomain.z0+k);
                            bufferIndex[index]=true;
                            numpoint++;
                        }
                    }
                }
            }
            for(int i=0;i<interpolation*Ny*Nz;i++){
                if(bufferIndex[i]){
                    fout<<(buffer[i].x-envelope)<<"\t"<<(buffer[i].y-envelope)<<"\t"<<(buffer[i].z-envelope)<<endl;
                    resultdata.get(buffer[i].x,buffer[i].y,buffer[i].z)=std::max(Nx,std::max(Ny,Nz));
                }
            }
            allcenter+=numpoint;
        }
        {
            int npx0 = datadomain.x0+(part-1)*interpolation;
            int npx1 = datadomain.x1;
            interpolation=npx1-npx0+1;
            std::vector<Dot3D> buffer(interpolation*Ny*Nz,Dot3D(0,0,0));
            std::vector<bool> bufferIndex(interpolation*Ny*Nz,false);
            int numpoint=0;
            #pragma omp parallel for reduction(+:numpoint)
            for(int i=0;i<interpolation;i++){
                for(int j=0;j<Ny;j++){
                    for(int k=0;k<Nz;k++){
                        int index=i*Ny*Nz+j*Nz+k;
                        if(maxpoint(Dot3D(npx0+i,datadomain.y0+j,datadomain.z0+k),datadomain,compareLength)){
                            buffer[index]=Dot3D(npx0+i,datadomain.y0+j,datadomain.z0+k);
                            bufferIndex[index]=true;
                            numpoint++;
                        }
                    }
                }
            }
            for(int i=0;i<interpolation*Ny*Nz;i++){
                if(bufferIndex[i]){
                    fout<<(buffer[i].x-envelope)<<"\t"<<(buffer[i].y-envelope)<<"\t"<<(buffer[i].z-envelope)<<endl;
                    resultdata.get(buffer[i].x,buffer[i].y,buffer[i].z)=std::max(Nx,std::max(Ny,Nz));
                }
            }
            allcenter+=numpoint;
        }
        fout.close();
        return allcenter;
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

    double porousDistribution(int n, int length){
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        RandomNumber rn;
        std::vector<double> buffer(n,0);
        DotList3D samples;
        while(samples.getN()<n){
            int tempx=rn.getRandomInt(datadomain.x0,datadomain.x1);
            int tempy=rn.getRandomInt(datadomain.y0,datadomain.y1);
            int tempz=rn.getRandomInt(datadomain.z0,datadomain.z1);

            if(contained(tempx+length-1, tempy+length-1, tempz+length-1, datadomain)){
                samples.addDot(Dot3D(tempx,tempy,tempz));
            }
        }

        double sum=0;
        #pragma omp parallel for reduction(+:sum)
        for(int i=0; i<n; i++){
            double oneans=calculateVoidage(samples.getDot(i), length);
            buffer[i]=oneans;
            sum+=oneans;
        }
        ofstream fout;
        fout.open("lvdata.dat", ios::ate);
        for(int i=0; i<n; i++){
            fout<<buffer[i]<<endl;
        }
        fout.close();
        return sum/n;
    }

    double smallWeightPath(Filter const& filter, addScale const& addscale, int sampleN){
        int samplenumber=sampleN;
        
        distanceTransform(filter, addscale ,2);
        Box3D datadomain=analysisdomain.enlarge(-envelope);
        
        //the transfer distance data become the rawdata
        rawdata=resultdata;

        resultdata.reset(1000000000);

        //find the max value of distance
        double maxvalue=0;
        for (int i = datadomain.x0; i <= datadomain.x1; i++){
            for (int j = datadomain.y0; j <= datadomain.y1; j++){
                for (int k = datadomain.z0; k <= datadomain.z1; k++){
                    maxvalue=std::max(maxvalue,rawdata.get(i,j,k));
                }
            }
        }

        DotList3D headsample;
        Box3D headdomain=datadomain.enlarge(-std::min(datadomain.getNx(),datadomain.getNy())/5);
        RandomNumber rn;
        while(headsample.getN()<samplenumber){
            int _x=rn.getRandomInt(headdomain.x0,headdomain.x1);
            int _y=rn.getRandomInt(headdomain.y0,headdomain.y1);
            int _z=datadomain.z0;
            if(abs(rawdata.get(_x,_y,_z)-0.0)>0.0001&&rawdata.get(_x,_y,_z)<maxvalue/3){
                headsample.addDot(Dot3D(_x,_y,_z));
            }
        }
        
        vector<std::stringstream> outputstream(samplenumber);
        vector<double> threshold(samplenumber,8*datadomain.getNz()*maxvalue);
        vector<DotList3D> pathans(samplenumber);
        vector<DotList3D> path(samplenumber);
        vector<double> pathvalueans(samplenumber,0);
        double pathweight=0;
        double inputvalue=1;

        omp_set_num_threads(samplenumber);
        #pragma omp parallel for
        for(int i=0;i<samplenumber;i++){
            ScalarField3D booldata=resultdata;
            Dot3D head=headsample.getDot(i);
            path[i].addDot(head);
            outputstream[i]<<"The head dot:"<<head.x<<","<<head.y<<","<<head.z<<std::endl;
            outputstream[i]<<"The head dot weight:"<<rawdata.get(head.x,head.y,head.z)<<std::endl;
            outputstream[i]<<"The biggest pore radius:"<<maxvalue<<std::endl;
            smallWeightDFS(pathweight, threshold[i], path[i], pathans[i], inputvalue, pathvalueans[i], Dot3D(head.x,head.y,head.z), maxvalue, addscale, filter, head, maxvalue, booldata);
            outputstream[i]<<"The path weight:"<<threshold[i]<<std::endl;
            outputstream[i]<<"The number of path dot:"<<pathans[i].getN()<<std::endl;
            outputstream[i]<<"The length of path:"<<pathvalueans[i]<<std::endl;
        }
        double sum=0;
        resultdata=rawdata;
        ofstream pathoutput;
        for(int i=0;i<samplenumber;i++){
            string name=createFileName("path",i,2);
            name+=".dat";
            pathoutput.open(name, ios::ate); 
            for(int j=0;j<pathans[i].getN();j++){
                Dot3D temp=pathans[i].getDot(j);
                resultdata.get(temp.x,temp.y,temp.z)=maxvalue*5;
                pathoutput<<temp.x<<"\t"<<temp.y<<"\t"<<temp.z<<endl;
            }
            pathoutput.close();
            sum+=pathvalueans[i];
            std::cout<<outputstream[i].str()<<endl;
        }
        return sum/samplenumber;
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
        
        //the change point
        addScale as(mode::d26);

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

    double calculateVoidage(Dot3D const& sample, int length){
        int sum=0;
        int total=pow(length,3);
        for(int i=sample.x; i<sample.x+length; i++){
            for(int j=sample.y; j<sample.y+length; j++){
                for(int k=sample.z; k<sample.z+length; k++){
                    int temp=rawdata.get(i,j,k);
                    sum+=temp;
                }
            }
        }

        return 1-static_cast<double>(sum)/total;
    }

    bool maxpoint(Dot3D const& center,Box3D const& datadomain,int compareLength){
        int centerx=center.x;
        int centery=center.y;
        int centerz=center.z;
        double centervalue=resultdata.get(centerx,centery,centerz);
        for(int i=-compareLength;i<=compareLength;i++){
            for(int j=-compareLength;j<=compareLength;j++){
                for(int k=-compareLength;k<=compareLength;k++){
                    if(!contained(centerx+i,centery+j,centerz+k,datadomain)) return false;
                    if(i==0&&j==0&&k==0) continue;
                    if(rawdata.get(centerx+i,centery+j,centerz+k)>0.001) return false;
                    if(resultdata.get(centerx+i,centery+j,centerz+k)>=centervalue) return false;
                }
            }
        }
        return true;
    }
    
    void smallWeightDFS(double pathweight, double& threshold, DotList3D& path, DotList3D& pathans, double pathvalue, double& pathvalueans,
                        Dot3D cur, double const& maxvalue, addScale const& addscale, Filter const& filter, Dot3D const& head, int domainlength, ScalarField3D & booldata){
        if(cur.z==analysisdomain.z1-envelope){
            if(pathweight<threshold){
                threshold=pathweight;
                pathans=path;
                pathvalueans=pathvalue;
                //printf("%f\n",threshold);
            }
            return;
        }
        if(!contained(cur.x,cur.y,cur.z,analysisdomain.enlarge(-envelope))) return;
        int radiuslength=std::pow(cur.x-head.x,2)+std::pow(cur.y-head.y,2);
        if(radiuslength>std::pow(domainlength,2)) return;
        if(pathweight>threshold) return;
        if(path.getN()>analysisdomain.getNz()+20) return;
        double curvalue=10*(maxvalue-rawdata.get(cur.x,cur.y,cur.z));
        double newvalue=curvalue+pathweight;
        //printf("%d,%d,%d,%f\n",cur.x,cur.y,cur.z,newvalue);
        if(newvalue>=booldata.get(cur.x,cur.y,cur.z)) return;
        booldata.get(cur.x,cur.y,cur.z)=newvalue;
        for(int i=0;i<addscale.n;i++){
            int addx=addscale.X[i],addy=addscale.Y[i],addz=addscale.Z[i];
            int x=cur.x+addx,y=cur.y+addy,z=cur.z+addz;
            double norm=filter.get(addx,addy,addz);
            path.addDot(Dot3D(cur.x,cur.y,cur.z));
            smallWeightDFS(newvalue,threshold, path, pathans, pathvalue+norm, pathvalueans,
                           Dot3D(x,y,z),maxvalue,addscale,filter,head,domainlength,booldata);
            path.eraseDot();
        }
    }
private:
    int envelope;
    Box3D analysisdomain;
    ScalarField3D rawdata;
    ScalarField3D resultdata;
};

#endif //ANALYSISTOOLS_H