#include "../cse452.h"
#include "ShapesUI.h"
#include "ShapesInterface.h"
#include "../Color.h"
#include <FL/Fl.H>
#include <FL/gl.h>
#include <GL/glu.h>
#include <iostream>
using namespace std;

ShapesUI::ShapesUI() {
    width = height = 0;

    // ToDo: initialize your variables here
}

ShapesUI::~ShapesUI() {
    // ToDo: delete your variables here
}

void ShapesUI::resize(int w, int h) {
    width = w;
    height = h;
}

void ShapesUI::draw() {
	// Sets up the viewport and background color
	setup3DDrawing(Color(0, 0, 0), width, height, true);

	// Changes the way triangles are drawn
	switch (shapesUI->getDisplayType()) {
	case DISPLAY_WIREFRAME: {
								glDisable(GL_LIGHTING);
								glPolygonMode(GL_FRONT, GL_LINE);
								glColor3f(1.0f, 1.0f, 1.0f);
	} break;
	case DISPLAY_FLAT_SHADING: {
								   glEnable(GL_LIGHTING);
								   glPolygonMode(GL_FRONT, GL_FILL);
								   glColor3f(1.0f, 1.0f, 1.0f);
								   glShadeModel(GL_FLAT);
	} break;
	case DISPLAY_SMOOTH_SHADING: {
									 glEnable(GL_LIGHTING);
									 glPolygonMode(GL_FRONT, GL_FILL);
									 glColor3f(1.0f, 1.0f, 1.0f);
									 glShadeModel(GL_SMOOTH);
	} break;
	default: break;
	}

	// Setup the camera
	gluLookAt(3.5 * cos(shapesUI->getYRot()) * cos(shapesUI->getXRot()),
		3.5 * sin(shapesUI->getYRot()),
		3.5 * cos(shapesUI->getYRot()) * sin(shapesUI->getXRot()), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	// ToDo: draw your shape here
	// DO NOT put the actual draw OpenGL code here - put it in the shape class and call the draw method
	
	shape->create();
	
	endDrawing();
}

int ShapesUI::handle(int event) {
    return 0;
}

void ShapesUI::changedShape()
{
	int tess1 = shapesUI->getTessel1();
	int tess2 = shapesUI->getTessel2();
	switch (shapesUI->getShapeType()) {
	case SHAPE_CUBE:
		shape = make_shared<Cube>(tess1);
		shape->create();
		break;
	case SHAPE_CYLINDER:
		shape = make_shared<Cylinder>(tess1, tess2);
		shape->create();
		break;
	case SHAPE_CONE:
		shape = make_shared<Cone>(tess1, tess2);
		shape->create();
		break;
	case SHAPE_SPHERE:
		shape = make_shared<Sphere>(tess1);
		shape->create();
		break;
	case SHAPE_HOURGLASS:
		shape = make_shared<Hourglass>(tess1, tess2);
		shape->create();
		break;
	}

    
    RedrawWindow();
}

void ShapesUI::changedTessel( ) {
   //essentially the same as ChangedShape()
	int tess1 = shapesUI->getTessel1();
	int tess2 = shapesUI->getTessel2();
	switch (shapesUI->getShapeType()) {
	case SHAPE_CUBE:
		shape = make_shared<Cube>(tess1);
		shape->create();
		break;
	case SHAPE_CYLINDER:
		shape = make_shared<Cylinder>(tess1, tess2);
		shape->create();
		break;
	case SHAPE_CONE:
		shape = make_shared<Cone>(tess1, tess2);
		shape->create();
		break;
	case SHAPE_SPHERE:
		shape = make_shared<Sphere>(tess1);
		shape->create();
		break;
	case SHAPE_HOURGLASS:
		shape = make_shared<Hourglass>(tess1, tess2);
		shape->create();
		break;
	}
    
    RedrawWindow();
}

