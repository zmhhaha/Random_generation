#ifndef TEMPLATE_H
#define TEMPLATE_H

//The template of all kind of point.
template<typename T>
class PointTemplate3D{
public:
    PointTemplate3D(T x_, T y_, T z_) : x(x_), y(y_), z(z_) { }
    PointTemplate3D():x(0),y(0),z(0){ }
    virtual ~PointTemplate3D(){ }
    void setPoint(T _x, T _y, T _z){
		x = _x;
		y = _y;
		z = _z;
	}
	void movePoint(T m_x, T m_y, T m_z) {
		x += m_x;
		y += m_y;
		z += m_z;
	};
    PointTemplate3D& operator+=(PointTemplate3D const& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
    PointTemplate3D& operator-=(PointTemplate3D const& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
public:
	T x, y, z;
};

#endif //TEMPLATE_H