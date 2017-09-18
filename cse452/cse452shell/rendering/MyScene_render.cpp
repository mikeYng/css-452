#include "../cse452.h"
#include "../sceneview/MyScene.h"
#include "RenderingInterface.h"
#include <FL/gl.h>
#include <cfloat>
using namespace std;

void MyScene::render(int type, int width, int height, unsigned char* pixels) {
    if (!isLoaded) {
        return;
    }

    // Add your rendering code here.
    // Keep track of your progress as a value between 0 and 1
    // so the progress bar can update as the rendering progresses
	rendering = true;
	progress = 0.0;
	int RECURSIVE_LIMIT = 5;
	Vector3 ray;
	Point3 pw; //point in the world
	Color col;
    switch (type) {
        case RenderingUI::RENDER_SCANLINE:  
			//dont need to do this
			progress = 1.0;
			break;
        case RenderingUI::RENDER_RAY_TRACING:  
			if (rendering) {
				for (int j = 0; j<height; j++){
					for (int i = 0; i<width; i++){
						pw = camera.getCameraToWorld() *  Point3((i + 0.5) * 2 / width - 1, (1 - (j + 0.5) * 2 / height), -1);
						ray = pw - camera.getEye();
						ray.normalize();
						col = calcColors(camera.getEye(), ray, RECURSIVE_LIMIT);
						//cout << col[0] << " " << col[1] << " " << col[2] << endl;
						putPixel(i, j, width, height, col, pixels);
					}
					progress = (j + 1)*1.0 / height;
					Fl::check();
					if (!rendering) {
						// break out of the loop if not rendering
						progress = 1.0;
						break;
					}
				}
			}
			
			break;
        case RenderingUI::RENDER_PATH_TRACING:  
			//dont need to do this
			progress = 1.0;
			break;
        default: break;
    }
}

void MyScene::stopRender()
{
	rendering = false;
	Fl::check();
    // Because this is threaded code, this function
    // can be called in the middle of your rendering code.
    // You should then stop at the next scanline
}

double MyScene::getRenderProgress() {
    // return the current progress as a value between 0 and 1
    return progress;
}

// add extra methods here

//essentially the same putPixel from the Brush lab
void MyScene::putPixel(int x, int y, int width, int height, Color & col, unsigned char * pixels) {
	int i = ((height - y - 1)*width + x) * 3;
	pixels[i] = (unsigned char)(col[0] * 255.0f);
	pixels[i + 1] = (unsigned char)(col[1] * 255.0f);
	pixels[i + 2] = (unsigned char)(col[2] * 255.0f);
}

void MyScene::rayIntersect(Point3 p, Vector3 dis, Point3 & hitP, Vector3 & hitN, Object * & obj) {
	double tmin = INT_MAX*1.0;
	double t, u, v;
	Point3 newHit;
	Vector3 newN;
	HitRecord hitObj;
	Matrix4 objtrans, invtrans, hittrans, invhittrans;

	for (vector<Object>::iterator it = objList.begin(); it != objList.end(); ++it){
		// I have no idea why I must use an iterator. It just wont work otherwise.
		//Object * current = it;
		objtrans = it->matrix;
		invtrans = objtrans.inverse();

		switch (it->shape_type)
		{
		case Object::CUBE:
			it->shape = make_shared<Cube>(5);
			hitObj = it->shape->intersect(invtrans*p, invtrans*dis);
			hitObj.sortHits();
			break;
		case Object::CYLINDER:
			it->shape = make_shared<Cylinder>(20,20);
			hitObj = it->shape->intersect(invtrans*p, invtrans*dis);
			hitObj.sortHits();
			break;
		case Object::CONE:
			it->shape = make_shared<Cone>(20, 20);
			hitObj = it->shape->intersect(invtrans*p, invtrans*dis);
			hitObj.sortHits();
			break;
		case Object::SPHERE:
			it->shape = make_shared<Sphere>(5);
			hitObj = it->shape->intersect(invtrans*p, invtrans*dis);
			hitObj.sortHits();
			break;
		}
		//first hit
		if (hitObj.getFirstHit(t, u, v, newHit, newN)){
			if (t>0 && t<tmin){
				tmin = t;
				hitP = newHit;
				hitN = newN;
				obj = &(*it);
				//cout << "in parsing"<<obj->diffuse[0] << endl;
				hittrans = objtrans;
				invhittrans = invtrans;
			}
		}

	}
	hitP = hittrans * hitP;
	hitN = hittrans * hitN;
	hitN.normalize();
}

Color MyScene::calcDiff(Object* obj, Point3 hitP, Vector3 hitN, Light* light) {
	Vector3 toLight = light->getPos() - hitP;
	bool inShadow = false;
	Object * current = 0;
	Point3 newhitP;
	Vector3 newhitN;
	rayIntersect(hitP + 0.001*unit(toLight), unit(toLight), newhitP, newhitN, current);
	if (current!=0) {
		if ((newhitP - hitP).length() < toLight.length())
			inShadow = true;
	}

	if (inShadow == false){
		Point3 fallOff = light->getFalloff();
		double attenuate = 1.0 / (fallOff[0] + fallOff[1] * toLight.length() + fallOff[2] * pow(toLight.length(), 2));
		double hitLight = hitN * unit(toLight);
		Color diffuse = obj->diffuse;
		//cout << "in diffuse calc " << diffuse[0] << endl;
		if (hitLight < 0) {
			hitLight = 0.0;
		}
		return (light->getColor())*diffuse*hitLight*attenuate;
	}
	else {
		return Color(0,0,0);
	}
	
}

Color MyScene::calcSpec(Object* obj, Vector3 ray, Point3 hitP, Vector3 hitN, Light* light) {
	Vector3 toLight = light->getPos() - hitP;
	bool inShadow = false;
	Object * current = 0;
	Point3 newhitP;
	Vector3 newhitN;
	rayIntersect(hitP + 0.001*unit(toLight), unit(toLight), newhitP, newhitN, current);
	if (current!=0) {
		if ((newhitP - hitP).length() < toLight.length())
			inShadow = true;
	}
	if (camera.getLook() * newhitN >0) {
		inShadow = true;
	}

	if (inShadow == false){
		Color lcol = light->getColor();
		Vector3 rv; //reflected vector
		rv = 2 * (toLight*hitN) * hitN - toLight;
		rv.normalize();
		Point3 fallOff = light->getFalloff();
		double attenuate = 1.0 / (fallOff[0] + fallOff[1] * toLight.length() + fallOff[2] * toLight.length());
		double rayReflect = ray * rv;
		double shine = obj->shine;
		// for some value, the shine can be out of these ranges, which is weird.
		if (shine <0 || shine >1) {
			shine = 1.0;
		}
		Color specular = obj->specular;
		//cout << "in spec calc" <<specular[0] << endl;
		if (rayReflect<0)
			rayReflect = 0.0;
		return lcol*specular*attenuate*pow(rayReflect, shine);
	}
	else {
		return Color(0,0,0);
	}
	
}

Color MyScene::calcColors(Point3 p, Vector3 ray, int level) {
	Color ambient = Color(0,0,0);
	Color diffuse = Color(0, 0, 0);
	Color reflect = Color(0, 0, 0);
	Color specular = Color(0, 0, 0);
	bool renderObj = true;

	Object * obj = 0; //set a new obj with a flag
	Point3 hitP;
	Vector3 hitN;
	rayIntersect(p, ray, hitP, hitN, obj);

	if (hitN*ray > 0) {
		renderObj = false;
	}
		

	if (obj!=0 && renderObj){

		ambient = obj->ambient;//trivial for ambient
        //calculate the diffuse and specular
		for (vector<Light>::iterator it = lights.begin(); it != lights.end(); ++it){
			diffuse = diffuse + calcDiff(obj, hitP, hitN, &(*it));
			specular = specular + calcSpec(obj, ray, hitP, hitN, &(*it));
		}

		//calculate reflection color
		Color oldReflect = obj->reflect;

		if (level>0 && oldReflect.getMax()>0){
			Vector3 newRay = ray;
			newRay = 2 * (newRay*hitN)*hitN - newRay;
			newRay.normalize();
			Color realRef = calcColors(hitP + 0.001*newRay, newRay, level - 1);
			reflect = realRef*oldReflect;
		}

		//final color
		//cout << ambient[0] << " " << specular[0] << " " << diffuse[0] << " " << reflect[0] << endl;
		Color finalCol = ambient + specular + diffuse + reflect;
		//normalize color
		double max = 1.0;
		for (int i = 0; i<3; i++){
			if (finalCol[i]<0) {
				finalCol[i] = 0.0;
			}
			else{
				if (finalCol[i]>max) {
					max = finalCol[i];
				}
			}
		}

		for (int j = 0; j<3; j++){
			if (max>1){
				finalCol[j] = finalCol[j] / max;
			}
		}
		return finalCol;
	}
	else {
		return background;
	}
		
}