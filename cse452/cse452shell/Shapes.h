#ifndef _SHAPES_H_
#define _SHAPES_H_
#include "../UIInterface.h"
#include "../cse452.h"
#include <FL/Fl_Window.H>
#include "../intersection/HitRecord.h"
#include "Point3.h"
#include "Vector3.h"
using namespace std;

class Shape {
public:
	Shape();
	void planeIntersect(Point3, Point3, Vector3, Vector3, HitRecord&);
	virtual void create()=0;
	virtual HitRecord intersect(Point3, Vector3) =0;
};


class Cube :public Shape {
public:
	int tess1;
	GLuint index;
	Cube(int);
	void create();
	HitRecord intersect(Point3, Vector3);
};

class Cone : public Shape {
public:
	int tess1;
	int tess2;
	GLuint index;
	Cone(int, int);
	void create();
	HitRecord intersect(Point3, Vector3);
};

class Cylinder : public Shape {
public:
	int tess1;
	int tess2;
	GLuint index;
	Cylinder(int, int);
	void create();
	HitRecord intersect(Point3, Vector3);
};

class Hourglass : public Shape {
public:
	int tess1;
	int tess2;
	GLuint index;
	Hourglass(int, int);
	void create();
	HitRecord intersect(Point3, Vector3);
	void drawFrame();
};

class Sphere : public Shape {
public:
	int tess1;
	GLuint index;
	Sphere(int );
	void create();
	HitRecord intersect(Point3, Vector3);
private:
	void subdivide(double*, double*, double*, int tess);
	void drawtri(double*, double*, double*);

};
#endif