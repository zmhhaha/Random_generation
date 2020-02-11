#ifndef FIBRESET_H
#define FIBRESET_H

#include<map>
#include"geometry.h"

class Fibre {
public:
    Fibre(){ }
    Fibre(Point3D const& point_,Vec3D const& vec_,double radius_) : point(point_),vec(vec_),radius(radius_) { }
    virtual ~Fibre(){}
    void setDirection(Vec3D const& vec_){
        vec=vec_;
    }
    void setBasePoint(Point3D const& point_){
        point=point_;
    }
    void setRadius(double radius_){
        radius=radius_;
    }
    UnitVec& getDirection(){
        return vec;
    }
    UnitVec const& getDirection() const{
        return vec;
    }
    Point3D& getBasePoint(){
        return point;
    }
    Point3D const& getBasePoint() const{
        return point;
    }
    double& getRadius(){
        return radius;
    }
    double const& getRadius() const{
        return radius;
    }
private:
    Point3D point;
    UnitVec vec;
    double radius;
};

class BetaFibre:public Fibre {
public:
    BetaFibre(Box3D const& box, map<double,double> const& rhomap, int n, double radius_, RandomNumber& rng){
        RandomPoint rp;
        Point3D basepoint=rp.NormalRandomPoint(box.x0,box.x1,box.y0,box.y1,box.z0,box.z1,rng);
        setBasePoint(basepoint);
		Point3D directionpoint=rp.BetaPoint(1000,rhomap,rng);
        UnitVec fibredirection(directionpoint);
		setDirection(fibredirection);
        setRadius(radius_);
    }
    virtual ~BetaFibre(){}
};

inline bool contained(int x, int y, int z, Fibre const& fibre) {
	Point3D point_origin(fibre.getBasePoint());
    Point3D spatial_point(x,y,z);
	Vec3D vec_origin(fibre.getDirection());
    Vec3D vec_two_point(point_origin,spatial_point);
    double radius=fibre.getRadius();
	double dot_vector = vec_origin.getNx()*vec_two_point.getNx()
                      + vec_origin.getNy()*vec_two_point.getNy()
                      + vec_origin.getNz()*vec_two_point.getNz();

	double X_shadow = point_origin.x + dot_vector*vec_origin.getNx();
	double Y_shadow = point_origin.y + dot_vector*vec_origin.getNy();
	double Z_shadow = point_origin.z + dot_vector*vec_origin.getNz();
	Point3D point_shadow(X_shadow, Y_shadow, Z_shadow);
	double distance = sqrt(pow(spatial_point.x - point_shadow.x, 2)
                         + pow(spatial_point.y - point_shadow.y, 2)
                         + pow(spatial_point.z - point_shadow.z, 2));
	if (distance <= radius) return true;
	else return false;
}

class FibreSet{
public:
    FibreSet(){
        fibreset.clear();
    }
    virtual ~FibreSet(){}
    void insertFibre(Fibre const& newfibre){
        fibreset.insert(std::pair<int,Fibre>(fibreset.size()+1,newfibre));
    }
    void earseFibre(int id){
        fibreset.erase(id);
    }
    std::map<int,Fibre>& getRawFibreSet(){
        return fibreset;
    }
    std::map<int,Fibre> const& getRawFibreSet() const{
        return fibreset;
    }
    Fibre & operator[](int id) {
        return fibreset.find(id)->second;
    }
    Fibre const& operator[](int id) const{
        return fibreset.find(id)->second;
    }
    int getN(){
        return fibreset.size();
    }
private:
    std::map<int,Fibre> fibreset; 
};

#endif //FIBRESET_H