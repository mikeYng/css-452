#include "../cse452.h"
#include "../Color.h"
#include <FL/Fl.H>
#include <FL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <Shapes.h>
#include "../Intersection/HitRecord.h"
using namespace std;

Shape::Shape() {
	//nothig 
}

void Shape::planeIntersect(Point3 x, Point3 p, Vector3 dir, Vector3 n, HitRecord & h) {
	double test = dot(dir, n);
	if (test != 0) {
		Vector3 xvec = Vector3(x[0], x[1], x[2]);
		Vector3 pvec = Vector3(p[0], p[1], p[2]);
		double t = (dot(xvec, n) - dot(n, pvec))/test;
		if (t > 0) {
			Point3 h0 = p + t*dir;
			if (n[0] != 0) {
				if ((abs(h0[1]) <= 0.5) && (abs(h0[2]) <= 0.5)) {
					// boundary checks, the intersection should be within the plane. 
					h.addHit(t, 0, 0, h0, n);
				}
			}
			else if (n[1] != 0) {
				if ((abs(h0[0]) <= 0.5) && (abs(h0[2]) <= 0.5)) {
					// boundary checks, the intersection should be within the plane. 
					h.addHit(t, 0, 0, h0, n);
				}
			}
			else {
				if ((abs(h0[0]) <= 0.5) && (abs(h0[1]) <= 0.5)) {
					// boundary checks, the intersection should be within the plane. 
					h.addHit(t, 0, 0, h0, n);
				}
			}
			
		}
	}
}
Cube::Cube(int tess) {
	this->tess1 = tess;
	this->index = glGenLists(1);
}

void Cube::create() {
	tess1 = max(tess1, 1); //cannot tessellate <1.
	double edge = 1.0 / tess1;
	glNewList(index, GL_COMPILE);
	glBegin(GL_TRIANGLES);
	// face #1
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(1, 0, 0);
			glVertex3f(0.5, 0.5 - j*edge, -0.5 + edge + i*edge);
			glVertex3f(0.5, 0.5 - edge - j*edge, -0.5 + i*edge);
			glVertex3f(0.5, 0.5 - j*edge, -0.5 + i*edge);
			glNormal3f(1, 0, 0);
			glVertex3f(0.5, 0.5 - j*edge, -0.5 + edge + i*edge);
			glVertex3f(0.5, 0.5 - edge - j*edge, -0.5 + edge + i*edge);
			glVertex3f(0.5, 0.5 - edge - j*edge, -0.5 + i*edge);
		}
	}
	//face #2 (back of #1)
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(-1, 0, 0);
			glVertex3f(-0.5, 0.5 - j*edge, 0.5 - edge - i*edge);
			glVertex3f(-0.5, 0.5 - edge - j*edge, 0.5 - i*edge);
			glVertex3f(-0.5, 0.5 - j*edge, 0.5 - i*edge);
			glNormal3f(-1, 0, 0);
			glVertex3f(-0.5, 0.5 - j*edge, 0.5 - edge - i*edge);
			glVertex3f(-0.5, 0.5 - edge - j*edge, 0.5 - edge - i*edge);
			glVertex3f(-0.5, 0.5 - edge - j*edge, 0.5 - i*edge);
		}
	}
	//face #3 (left of face 1)
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(0, 0, 1);
			glVertex3f(0.5 - edge - i*edge, 0.5 - j*edge, 0.5);
			glVertex3f(0.5 - i*edge, 0.5 - edge - j*edge, 0.5);
			glVertex3f(0.5 - i*edge, 0.5 - j*edge, 0.5);
			glNormal3f(0, 0, 1);
			glVertex3f(0.5 - edge - i*edge, 0.5 - j*edge, 0.5);
			glVertex3f(0.5 - edge - i*edge, 0.5 - edge - j*edge, 0.5);
			glVertex3f(0.5 - i*edge, 0.5 - edge - j*edge, 0.5);
		}
	}
	//face #4 (right of #1, back of #3)
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(0, 0, -1);
			glVertex3f(-0.5 + edge + i*edge, 0.5 - j*edge, -0.5);
			glVertex3f(-0.5 + i*edge, 0.5 - edge - j*edge, -0.5);
			glVertex3f(-0.5 + i*edge, 0.5 - j*edge, -0.5);
			glNormal3f(0, 0, -1);
			glVertex3f(-0.5 + edge + i*edge, 0.5 - j*edge, -0.5);
			glVertex3f(-0.5 + edge + i*edge, 0.5 - edge - j*edge, -0.5);
			glVertex3f(-0.5 + i*edge, 0.5 - edge - j*edge, -0.5);
		}
	}
	//face #5(upper face)
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(0, 1, 0);
			glVertex3f(-0.5 + i*edge, 0.5, -0.5 + edge + j*edge);
			glVertex3f(-0.5 + i*edge + edge, 0.5, -0.5 + j*edge);
			glVertex3f(-0.5 + i*edge, 0.5, -0.5 + j*edge);
			glNormal3f(0, 1, 0);
			glVertex3f(-0.5 + i*edge, 0.5, -0.5 + edge + j*edge);
			glVertex3f(-0.5 + edge + i*edge, 0.5, -0.5 + edge + j*edge);
			glVertex3f(-0.5 + edge + i*edge, 0.5, -0.5 + j*edge);
		}
	}
	//face #6(lower face)
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess1; ++j) {
			glNormal3f(0, -1, 0);
			glVertex3f(0.5 - i*edge, -0.5, -0.5 + edge + j*edge);
			glVertex3f(0.5 - i*edge - edge, -0.5, -0.5 + j*edge);
			glVertex3f(0.5 - i*edge, -0.5, -0.5 + j*edge);
			glNormal3f(0, -1, 0);
			glVertex3f(0.5 - i*edge, -0.5, -0.5 + edge + j*edge);
			glVertex3f(0.5 - edge - i*edge, -0.5, -0.5 + edge + j*edge);
			glVertex3f(0.5 - i*edge - edge, -0.5, -0.5 + j*edge);
		}
	}
	glEnd();
	glEndList();
	glCallList(index);
	glDeleteLists(index, 1);
}

HitRecord Cube::intersect(Point3 p, Vector3 dir) {
	HitRecord h;
	planeIntersect(Point3(0.5, 0, 0), p, dir, Vector3(1, 0, 0), h);
	planeIntersect(Point3(-0.5, 0, 0), p, dir, Vector3(-1, 0, 0),h);
	planeIntersect(Point3(0, 0.5, 0), p, dir, Vector3(0, 1, 0),h);
	planeIntersect(Point3(0, -0.5, 0), p, dir, Vector3(0, -1, 0),h);
	planeIntersect(Point3(0, 0, 0.5), p, dir, Vector3(0, 0, 1),h);
	planeIntersect(Point3(0, 0, -0.5), p, dir, Vector3(0, 0, -1),h);
	return h;
}
Cylinder::Cylinder(int tess1, int tess2) {
	this->tess1 = tess1;
	this->tess2 = tess2;
	this->index = glGenLists(1);
}

void Cylinder::create() {
	tess1 = max(tess1, 3); //cannot tessellate <3.
	tess2 = max(tess2, 1); //cannot tessellate <1.
	double edge2 = 1.0 / tess2;
	glNewList(index, GL_COMPILE);
	glBegin(GL_TRIANGLES);
	//top surface
	for (int i = 0; i < tess1; ++i) {
		glNormal3f(0, 1, 0);
		glVertex3f(0.5*cos(2 * i*M_PI / tess1), 0.5, 0.5*sin(2 * i*M_PI / tess1));
		glVertex3f(0, 0.5, 0);
		glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0.5, 0.5*sin(2 * (i + 1)*M_PI / tess1));
	}
	//bottom surface
	for (int i = 0; i < tess1; ++i) {
		glNormal3f(0, -1, 0);
		glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), -0.5, 0.5*sin(2 * (i + 1)*M_PI / tess1));
		glVertex3f(0, -0.5, 0);
		glVertex3f(0.5*cos(2 * i*M_PI / tess1), -0.5, 0.5*sin(2 * i*M_PI / tess1));
	}
	//the side surfaces
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess2; ++j) {
			glNormal3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0, 0.5*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0.5 - j*edge2, 0.5*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0.5 - edge2 - j*edge2, 0.5*sin(2 * (i + 1)*M_PI / tess1));
			glNormal3f(0.5*cos(2 * i*M_PI / tess1), 0, 0.5*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5*cos(2 * i*M_PI / tess1), 0.5 - j*edge2, 0.5*sin(2 * i*M_PI / tess1));

			glNormal3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0, 0.5*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), 0.5 - edge2 - j*edge2, 0.5*sin(2 * (i + 1)*M_PI / tess1));
			glNormal3f(0.5*cos(2 * i*M_PI / tess1), 0, 0.5*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5*cos(2 * i*M_PI / tess1), 0.5 - edge2 - j*edge2, 0.5*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5*cos(2 * i*M_PI / tess1), 0.5 - j*edge2, 0.5*sin(2 * i*M_PI / tess1));
		}
	}
	glEnd();
	glEndList();
	//cout << "in shape: " << index << endl;
	glCallList(index);
	glDeleteLists(index, 1);
}

HitRecord Cylinder::intersect(Point3 p, Vector3 dir) {
	HitRecord h;
	// the side cylinder
	double a = pow(dir[0], 2) + pow(dir[2], 2);
	double b = 2.0*(p[0] * dir[0] + p[2] * dir[2]);
	double c = pow(p[0], 2) + pow(p[2], 2) - 0.25;
	double test = pow(b, 2) - 4.0*a*c;
	if (test < 0) {
		return h;
	}
	else if (test == 0) {
		double t = -b / (2.0*a);
		Point3 h0 = p + t*dir;
		Vector3 n = Vector3(h0[0], 0, h0[2]);
		n.normalize();
		h.addHit(t, 0, 0, h0, n);
	}
	else {
		double t1 = (-b - sqrt(test)) / (2.0*a);
		double t2 = (-b + sqrt(test)) / (2.0*a);
		if (t1 > 0) {
			Point3 h1 = p + t1*dir;
			if (abs(h1[1]) <= 0.5) {
				Vector3 n1 = Vector3(h1[0], 0, h1[2]);
				n1.normalize();
				h.addHit(t1, 0, 0, h1, n1);
			}
			
		}
		if (t2 > 0) {
			Point3 h2 = p + t2*dir;
			if (abs(h2[1]) <= 0.5) {
				Vector3 n2 = Vector3(h2[0], 0, h2[2]);
				n2.normalize();
				h.addHit(t2, 0, 0, h2, n2);
			}
			
		}
	}

	// the caps
	double topt = (0.5 + p[1]) / (-dir[1]);
	if (topt > 0) {
		Point3 htop = p + topt*dir;
		if ((pow(htop[0], 2) + pow(htop[2], 2)) <= 0.25) {
			h.addHit(topt, 0, 0, htop, Vector3(0, -1, 0));
		}
	}
	
	double bott = (-0.5 + p[1]) / (-dir[1]);
	if (bott > 0) {
		Point3 hbot = p + bott*dir;
		if ((pow(hbot[0], 2) + pow(hbot[2], 2)) <= 0.25) {
			h.addHit(bott, 0, 0, hbot, Vector3(0, 1, 0));
		}
	}
	
	return h;
}

Cone::Cone(int tess1, int tess2) {
	this->tess1 = tess1;
	this->tess2 = tess2;
	this->index = glGenLists(1);
}

void Cone::create() {
	glNewList(index, GL_COMPILE);
	glBegin(GL_TRIANGLES);
	tess1 = max(tess1, 3); //cannot tessellate <3.
	tess2 = max(tess2, 1); //cannot tessellate <1.
	//bottom surface
	for (int i = 0; i < tess1; ++i) {
		glNormal3f(0, -1, 0);
		glVertex3f(0.5*cos(2 * (i + 1)*M_PI / tess1), -0.5, 0.5*sin(2 * (i + 1)*M_PI / tess1));
		glVertex3f(0, -0.5, 0);
		glVertex3f(0.5*cos(2 * i*M_PI / tess1), -0.5, 0.5*sin(2 * i*M_PI / tess1));
	}
	//sides
	//the first one is special
	for (int i = 0; i < tess1; ++i)  {
		double x_i = 0.5 / tess2*cos(2 * i*M_PI / tess1);
		double z_i = 0.5 / tess2*sin(2 * i*M_PI / tess1);
		double x_i1 = 0.5 / tess2*cos(2 * (i + 1)*M_PI / tess1);
		double z_i1 = 0.5 / tess2*sin(2 * (i + 1)*M_PI / tess1);
		glNormal3f(0, 0, 0);
		glVertex3f(0, 0.5, 0);
		glNormal3f(x_i1, 0.5*sqrt(pow(x_i1, 2) + pow(z_i1, 2)), z_i1);
		glVertex3f(x_i1, 0.5 - 1.0 / tess2, z_i1);
		glNormal3f(x_i, 0.5*sqrt(pow(x_i, 2) + pow(z_i, 2)), z_i);
		glVertex3f(x_i, 0.5 - 1.0 / tess2, z_i);
	}
	//the other ones
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess2 - 1; ++j) {
			double ratio = (j + 1)*1.0 / tess2;
			double nex_ratio = (j + 2)*1.0 / tess2;
			glNormal3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5 - 1.0*ratio, 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
			glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5 - 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
			glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5 - 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));

			glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5 - 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));
			glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
			glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5 - 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
			glNormal3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
			glVertex3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5 - 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
		}
	}
	glEnd();
	glEndList();
	//cout << "in shape: " << index << endl;
	glCallList(index);
	glDeleteLists(index, 1);
}

HitRecord Cone::intersect(Point3 p, Vector3 dir) {
	HitRecord h;
	double a = pow(dir[0], 2) + pow(dir[2], 2) - 0.25*pow(dir[1], 2);
	double b = 2 * p[0] * dir[0] + 2 * p[2] * dir[2] + 0.5*dir[1] * (0.5 - p[1]);
	double c = pow(p[0], 2) + pow(p[2], 2) - 0.25*pow((0.5 - p[1]), 2);
	double test = pow(b, 2) - 4.0*a*c;
	if (test < 0) {
		return h;
	}
	else if (test == 0) {
		double t = -b / (2.0*a);
		Point3 h0 = p + t*dir;
		Vector3 n = Vector3(h0[0], 0, h0[2]);
		n.normalize();
		h.addHit(t, 0, 0, h0, n);
	}
	else {
		double t1 = (-b - sqrt(test)) / (2.0*a);
		double t2 = (-b + sqrt(test)) / (2.0*a);
		if (t1 > 0) {
			Point3 h1 = p + t1*dir;
			if (abs(h1[1]) <= 0.5) {
				Vector3 n1 = Vector3(h1[0], 0.5*sqrt(pow(h1[0],2) + pow(h1[2], 2)), h1[2]);
				n1.normalize();
				h.addHit(t1, 0, 0, h1, n1);
			}

		}
		if (t2 > 0) {
			Point3 h2 = p + t2*dir;
			if (abs(h2[1]) <= 0.5) {
				Vector3 n2 = Vector3(h2[0], 0.5*sqrt(pow(h2[0], 2) + pow(h2[2], 2)), h2[2]);
				n2.normalize();
				h.addHit(t2, 0, 0, h2, n2);
			}

		}
	}
	// the caps
	double topt = (0.5 + p[1]) / (-dir[1]);
	if (topt > 0) {
		Point3 htop = p + topt*dir;
		if ((pow(htop[0], 2) + pow(htop[2], 2)) <= 0.25) {
			h.addHit(topt, 0, 0, htop, Vector3(0, -1, 0));
		}
	}
	
	return h;
}
Hourglass::Hourglass(int tess1, int tess2) {
	this->tess1 = tess1;
	this->tess2 = tess2;
	this->index = glGenLists(1);
}


void Hourglass::drawFrame() {
	double length = 0.2;
	// top bar
	for (int i = 0; i < 5; ++i) {
		glNormal3f(0, 1, 0);
		glVertex3f(-length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, 0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, 0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, 0.6, -0.5 + i*length);
	}
	//need to do a second time to make sure it shades both sides. 
	for (int i = 0; i < 5; ++i) {
		glNormal3f(0, -1, 0);
		glVertex3f(length / 2.0, 0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(-length / 2.0, 0.6, -0.5 + i*length);
		glVertex3f(length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(-length / 2.0, 0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, 0.6, -0.5 + i*length);
	}

	// bot bar
	for (int i = 0; i < 5; ++i) {
		glNormal3f(0, 1, 0);
		glVertex3f(-length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, -0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, -0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, -0.6, -0.5 + i*length);
	}
	for (int i = 0; i < 5; ++i) {
		glNormal3f(0, -1, 0);
		glVertex3f(length / 2.0, -0.6, -0.5 + i*length);
		glVertex3f(-length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(-length / 2.0, -0.6, -0.5 + i*length);
		glVertex3f(length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(-length / 2.0, -0.6, -0.5 + length + i*length);
		glVertex3f(length / 2.0, -0.6, -0.5 + i*length);
	}

	//left side bar
	for (int i = 0; i < 6; ++i) {
		glNormal3f(0, 0, 1);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, 0.5);
		glVertex3f(-length / 2.0, 0.6 - i*length, 0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, 0.5);
	}
	for (int i = 0; i < 6; ++i) {
		glNormal3f(0, 0, -1);
		glVertex3f(length / 2.0, 0.6 - i*length, 0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(-length / 2.0, 0.6 - i*length, 0.5);
		glVertex3f(length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, 0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, 0.5);
	}

	//right side bar
	for (int i = 0; i < 6; ++i) {
		glNormal3f(0, 0, 1);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, -0.5);
		glVertex3f(-length / 2.0, 0.6 - i*length, -0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, -0.5);
	}
	for (int i = 0; i < 6; ++i) {
		glNormal3f(0, 0, -1);
		glVertex3f(length / 2.0, 0.6 - i*length, -0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(-length / 2.0, 0.6 - i*length, -0.5);
		glVertex3f(length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(-length / 2.0, 0.6 - length - i*length, -0.5);
		glVertex3f(length / 2.0, 0.6 - i*length, -0.5);
	}
}
void Hourglass::create() {
	glNewList(index, GL_COMPILE);
	glBegin(GL_TRIANGLES);
	tess1 = max(tess1, 3); //cannot tessellate <3.
	tess2 = max(tess2, 1); //cannot tessellate <1.
	//the frame to hold the hourglass
	drawFrame();
	//bottom surface, same as cone
	for (int i = 0; i < tess1; ++i) {
		glNormal3f(0, -1, 0);
		glVertex3f(0.25*cos(2 * (i + 1)*M_PI / tess1), -0.5, 0.25*sin(2 * (i + 1)*M_PI / tess1));
		glVertex3f(0, -0.5, 0);
		glVertex3f(0.25*cos(2 * i*M_PI / tess1), -0.5, 0.25*sin(2 * i*M_PI / tess1));
	}
	// top surface
	for (int i = 0; i < tess1; ++i) {
		glNormal3f(0, 1, 0);
		glVertex3f(0, 0.5, 0);
		glVertex3f(0.25*cos(2 * (i + 1)*M_PI / tess1), 0.5, 0.25*sin(2 * (i + 1)*M_PI / tess1));
		glVertex3f(0.25*cos(2 * i*M_PI / tess1), 0.5, 0.25*sin(2 * i*M_PI / tess1));
	}
	//sides
	//the first one is special
	for (int i = 0; i < tess1; ++i)  {
		double x_i = 0.25 / tess2*cos(2 * i*M_PI / tess1);
		double z_i = 0.25 / tess2*sin(2 * i*M_PI / tess1);
		double x_i1 = 0.25 / tess2*cos(2 * (i + 1)*M_PI / tess1);
		double z_i1 = 0.25 / tess2*sin(2 * (i + 1)*M_PI / tess1);
		glNormal3f(x_i1, 0.5*sqrt(pow(x_i1, 2) + pow(z_i1, 2)), z_i1);
		glVertex3f(x_i1, 0 - 0.5 / tess2, z_i1);
		glNormal3f(0, 0, 0);
		glVertex3f(0, 0, 0);
		glNormal3f(x_i, 0.5*sqrt(pow(x_i, 2) + pow(z_i, 2)), z_i);
		glVertex3f(x_i, 0 - 0.5 / tess2, z_i);
	}
	//need two of that
	for (int i = 0; i < tess1; ++i)  {
		double x_i = 0.25 / tess2*cos(2 * i*M_PI / tess1);
		double z_i = 0.25 / tess2*sin(2 * i*M_PI / tess1);
		double x_i1 = 0.25 / tess2*cos(2 * (i + 1)*M_PI / tess1);
		double z_i1 = 0.25 / tess2*sin(2 * (i + 1)*M_PI / tess1);
		glNormal3f(x_i1, 0.5*sqrt(pow(x_i1, 2) + pow(z_i1, 2)), z_i1);
		glVertex3f(x_i1, 0 + 0.5 / tess2, z_i1);
		glNormal3f(0, 0, 0);
		glVertex3f(0, 0, 0);
		glNormal3f(x_i, 0.5*sqrt(pow(x_i, 2) + pow(z_i, 2)), z_i);
		glVertex3f(x_i, 0 + 0.5 / tess2, z_i);
	}
	//the other ones, this is the bottom cone
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess2 - 1; ++j) {
			double ratio = (j + 1)*0.5 / tess2;
			double nex_ratio = (j + 2)*0.5 / tess2;
			if (j > tess2 / 2.0) {
				glNormal3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));

				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
			}
			else {
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));

				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), -1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), -1.0*nex_ratio, 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
			}

		}
	}
	//the top sides
	for (int i = 0; i < tess1; ++i) {
		for (int j = 0; j < tess2 - 1; ++j) {
			double ratio = (j + 1)*0.5 / tess2;
			double nex_ratio = (j + 2)*0.5 / tess2;
			if (j > tess2 / 2.0) {
				glNormal3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));

				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
			}
			else {
				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));

				glNormal3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * (i + 1)*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * (i + 1)*M_PI / tess1));
				glNormal3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *ratio*cos(2 * i*M_PI / tess1), 1.0*ratio, 0.5*ratio*sin(2 * i*M_PI / tess1));
				glNormal3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 0.5*sqrt(pow(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 2) + pow(0.5*nex_ratio*sin(2 * i*M_PI / tess1), 2)), 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
				glVertex3f(0.5 *nex_ratio*cos(2 * i*M_PI / tess1), 1.0*nex_ratio, 0.5*nex_ratio*sin(2 * i*M_PI / tess1));
			}
			
		}
	}
	glEnd();
	glEndList();
	glCallList(index);
	glDeleteLists(index, 1);
}

HitRecord Hourglass::intersect(Point3 p, Vector3 dir) {
	HitRecord h;
	return h;
}
Sphere::Sphere(int tess1) {
	this->tess1 = tess1;
	this->index = glGenLists(1);
}


void Sphere::drawtri(double* p1, double* p2, double* p3) {
	glBegin(GL_TRIANGLES);
	glNormal3f(p1[0], p1[1], p1[2]);
	glVertex3f(p1[0], p1[1], p1[2]);
	glNormal3f(p2[0], p2[1], p2[2]);
	glVertex3f(p2[0], p2[1], p2[2]);
	glNormal3f(p3[0], p3[1], p3[2]);
	glVertex3f(p3[0], p3[1], p3[2]);
	glEnd();
}


void Sphere::subdivide(double* p1, double* p2, double* p3, int tess) {
	if (tess < 1) {
		tess = 1;
	}
	if (tess == 1) {
		drawtri(p1, p2, p3);
	}
	else {
		double mids[3][3];
		double* points[4] = { p1, p2, p3, p1 };
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				mids[i][j] = (points[i][j] + points[i + 1][j]) / 2.0;
			}
		}
		for (int i = 0; i < 3; ++i) {
			double norm = sqrt(pow(mids[i][0], 2) + pow(mids[i][1], 2) + pow(mids[i][2], 2));
			for (int j = 0; j < 3; ++j) {
				mids[i][j] = mids[i][j] / norm *0.5;
			}
		}
		subdivide(p1, mids[0], mids[2], tess - 1);
		subdivide(mids[0], p2, mids[1], tess - 1);
		subdivide(mids[2], mids[1], p3, tess - 1);
		subdivide(mids[0], mids[1], mids[2], tess - 1);

	}
}

void Sphere::create() {
	//values specially calculated such that with tess=1, these coordinates would generate a distance of 0.5
	tess1 = min(tess1, 5);
	glNewList(index, GL_COMPILE);
	double val1 = 0.525731112119133606/2.0;
	double val2 = 0.850650808352039932/2.0;
	double p0[3] = { -val1, 0, val2 };
	double p1[3] = { val1, 0, val2 };
	double p2[3] = { -val1, 0, -val2 };
	double p3[3] = { val1, 0, -val2 };
	double p4[3] = { 0, val2, val1 };
	double p5[3] = { 0, val2, -val1 };
	double p6[3] = { 0, -val2, val1 };
	double p7[3] = { 0, -val2, -val1 };
	double p8[3] = { val2, val1, 0 };
	double p9[3] = { -val2, val1, 0 };
	double p10[3] = { val2, -val1, 0 };
	double p11[3] = { -val2, -val1, 0 };

	subdivide(p1, p4, p0, tess1);
	subdivide(p4, p9, p0, tess1);
	subdivide(p4, p5, p9, tess1);
	subdivide(p8, p5, p4, tess1);
	subdivide(p1, p8, p4, tess1);
	subdivide(p1, p10, p8, tess1);
	subdivide(p10, p3, p8, tess1);
	subdivide(p8, p3, p5, tess1);
	subdivide(p3, p2, p5, tess1);
	subdivide(p3, p7, p2, tess1);
	subdivide(p3, p10, p7, tess1);
	subdivide(p10, p6, p7, tess1);
	subdivide(p6, p11, p7, tess1);
	subdivide(p6, p0, p11, tess1);
	subdivide(p6, p1, p0, tess1);
	subdivide(p10, p1, p6, tess1);
	subdivide(p11, p0, p9, tess1);
	subdivide(p2, p11, p9, tess1);
	subdivide(p5, p2, p9, tess1);
	subdivide(p11, p2, p7, tess1);
	glEndList();
	glCallList(index);
	glDeleteLists(index, 1);
}

HitRecord Sphere::intersect(Point3 p, Vector3 dir) {
	HitRecord h;
	double a = pow(dir[0], 2) + pow(dir[1], 2) + pow(dir[2], 2);
	double b = 2.0 * (dir[0] * p[0] + dir[1] * p[1] + dir[2] * p[2]);
	double c = pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2) - 0.25;
	double test = pow(b, 2) - 4 * a*c;
	if (test < 0) {
		return h;
	}
	else if (test == 0) {
		double t = -b  / (2.0*a);
		if (t > 0) {
			Point3 h0 = p + t*dir;
			Vector3 n = Vector3(h0[0], h0[1], h0[2]);
			n.normalize();
			h.addHit(t, 0, 0, h0, n);
		}
	}
	else {
		double t1 = (-b - sqrt(test)) / (2.0*a);
		double t2 = (-b + sqrt(test)) / (2.0*a);
		if (t1 > 0) {
			Point3 h1 = p + t1*dir;
			Vector3 n1 = Vector3(h1[0], h1[1], h1[2]);
			n1.normalize();
			h.addHit(t1, 0, 0, h1, n1);
		}
		if (t2 > 0) {
			Point3 h2 = p + t2*dir;
			Vector3 n2 = Vector3(h2[0], h2[1], h2[2]);
			n2.normalize();
			h.addHit(t2, 0, 0, h2, n2);
		}
	}
	return h;
}