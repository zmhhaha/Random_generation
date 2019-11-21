#ifndef VEC_H
#define VEC_H

#include<cmath>
#include<vector>
#include"geometry.h"
#include"template.h"
#include"randomnumber.h"

typedef PointTemplate3D<double> Point3D;

class RandomPoint: public Point3D{
public:
    RandomPoint(){}
    virtual ~RandomPoint(){}
    Point3D& NormalRandomPoint(double downx,double upx,
                               double downy,double upy,
                               double downz,double upz,
                               RandomNumber& rng){
        x=rng.getRandomfloat(downx,upx);
        y=rng.getRandomfloat(downy,upy);
        z=rng.getRandomfloat(downz,upz);
        return *this;
    }
    Point3D& BetaPoint(int n, map<double,double> dis, RandomNumber& rng){
        x=rng.getRandomfloat(-n,n);
        y=rng.getRandomfloat(-n,n);
        z=tan(rng.getRandomfromDistribution(dis)-PI/2)*sqrt(pow(x,2)+pow(y,2));
        return *this;
    }
};

class Vec3D{
public:
    Vec3D(Point3D const& begin_, Point3D const& end_):begin(begin_),end(end_) { }
    Vec3D(Point3D const& end_):begin(Point3D()),end(end_) { }
    Vec3D(double x0_, double y0_,double z0_,
          double x1_, double y1_, double z1_)
          :begin(x0_, y0_, z0_),end(x1_, y1_, z1_) { }
    Vec3D():begin(Point3D()),end(Point3D(1,1,1)) { }
    virtual ~Vec3D() { }
    double norm() const{
        return sqrt(pow(end.x-begin.x,2)+pow(end.y-begin.y,2)+pow(end.z-begin.z,2));
    }
    double getNx(){return end.x-begin.x;}
    double getNy(){return end.y-begin.y;}
    double getNz(){return end.z-begin.z;}
public:
    Point3D begin,end;
};

class UnitVec:public Vec3D {
public:
    UnitVec(Vec3D const& vc_){
        iniVec(vc_.begin,vc_.end);
    }
    UnitVec(Point3D const& begin_, Point3D const& end_){
        iniVec(begin_,end_);
    }
    UnitVec(Point3D const& end_){
        iniVec(Point3D(0,0,0),end_);
    }
    UnitVec(){
        iniVec(Point3D(0,0,0),Point3D(1,1,1));
    }
    virtual ~UnitVec(){};
    void iniVec(Point3D const& begin_, Point3D const& end_){
        Vec3D vc(begin_,end_);
        begin=Point3D();
        end=Point3D( vc.getNx() /vc.norm(), vc.getNy() /vc.norm(), vc.getNz() /vc.norm());
    }
};



#endif //VEC_H