#include "../cse452.h"
#include "MyScene.h"

void MyScene::resize(int w, int h) {
    // resize the film plane to the specified width/height
    camera.setWidthHeight(w, h);
}

/// Note: your camera and screen clear, etc, will be set up by
/// SceneviewUI.cpp *before* this gets called
void MyScene::draw() {
    // render the scene using OpenGL
    if (!isLoaded) // Don't draw if loadSceneFile hasn't been called yet
        return;


    // Turn off all lights
    for ( int i = 0; i < 7; i++ )
        glDisable( GL_LIGHT0 + i );

    //  .. and reset
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, &ambientLight[0] );
    for (unsigned int i = 0; i < lights.size(); i++) {
        lights[i].setOpenGLLight( GL_LIGHT0 + i );
    }

    // TODO: draw the rest of the scene here
	for (unsigned i = 0; i < objList.size(); ++i) {
		Object current = objList[i];
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, &current.ambient[0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &current.diffuse[0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &current.specular[0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, &current.emitted[0]);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, current.shine);
		glPushMatrix();
		glMultMatrixd(&objList[i].matrix(0,0));
		switch (objList[i].shape_type) {
		case Object::CUBE:
			current.shape = make_shared<Cube>(5);
			current.shape->create();
			break;
		case Object::CYLINDER:
			current.shape = make_shared<Cylinder>(20, 20);
			current.shape->create();
			break;
		case Object::CONE:
			current.shape = make_shared<Cone>(20, 20);
			current.shape->create();
			break;
		case Object::SPHERE:
			current.shape = make_shared<Sphere>(5);
			current.shape->create();
			break;
		}
	
		glPopMatrix();
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glColor4f(current.transparent[0], current.transparent[1], current.transparent[2], 0.5);
		//delete objList[i];
	}
}
