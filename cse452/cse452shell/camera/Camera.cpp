#include "../cse452.h"
#include "Camera.h"
#include <cmath>
#include <FL/Fl.H>
using namespace std;

Camera::Camera() 
{
	p = Point3 (3, 2, 6);
	L = Vector3(-3, -2, -6);
	U = Vector3(0, 1, 0);
	n = unit(-L);
	v = unit(U - dot(U, n)*n);
	u = unit(cross(v, n)); 
	dn = 0.01;
	df = 1000000;
	height = 5;
	width = 5; 
	zoom = 0;
	setTranslation();
	setRotation();
	setWidthScaling();
	setFarScaling();
	setPerspective();

}

Camera::~Camera() {
    // destroy your data here
}

// The following three should be unit length and orthogonal to each other
// u vector
Vector3 Camera::getRight() const
{
    return u;
}

// v vector
Vector3 Camera::getUp() const
{
    return v;
}

// - n vector
Vector3 Camera::getLook() const
{
    return -n;
}

double Camera::getSkew() const
{
    // Change this to implement the extra credit
    return 0.0;
}

double Camera::getAspectRatioScale() const
{
    // Change this to implement the extra credit
    return 1.0;
}

Point3 Camera::getProjectionCenter() const
{
    // Change this to implement the extra credit
    return Point3( 0.0, 0.0, 0.0 );
}

Matrix4 Camera::getProjection() const {
    // return the current projection and scale matrix
	//cout << "projection" << D*Sxyz*Sxy << endl;
    return D*Sxyz*Sxy;
}

Matrix4 Camera::getWorldToCamera() const {
    // return the current world to camera matrix
    // Rotation and translation
	//cout << "world to camera" << R*T << endl;
    return R*T;
}

Matrix4 Camera::getRotationFromXYZ() const
{
    // return just the rotation matrix
    return R;
}

Matrix4 Camera::getRotationToXYZ() const
{
    // return just the rotation matrix

    return R_1;
}

Matrix4 Camera::getCameraToWorld() const {
    // return the current camera to world matrix
    // This is the inverse of the rotation, translation, and scale

    // Change this
    return T_1*R_1*Sxy_1*Sxyz_1;
}

int Camera::getWidth()  const{
    // return the current image width
    return width;
}

int Camera::getHeight()  const{
    // return the current image height
    return height;
}

Point3 Camera::getEye()  const{
    // return the current eye location
    return p;
}

double Camera::getZoom() const
{
    return zoom;
}

void Camera::setFrom(const Point3& from) {
	p = from;
	setTranslation();
}

void Camera::setAt(const Point3& at) {
	L = at - p;
	setRotation();
}

void Camera::setLook(const Vector3& l) {
	L = l;
	setRotation();
}

void Camera::setUp(const Vector3& up) {
	U = up;
	setRotation();
}

void Camera::setWidthHeight(int w, int h) {
	width = w;
	height = h;
	setWidthScaling();
}

void Camera::setZoom(double z) {
	zoom = z;
	setWidthScaling();
    
}

void Camera::setNearFar(double n, double f) {
	dn = n;
	df = f;
	setFarScaling();
	setPerspective();
}

void Camera::setSkew( double d )
{
}

void Camera::setAspectRatioScale( double d )
{
}

void Camera::setProjectionCenter( const Point3 &p )
{
}

void Camera::moveForward(double dist) {
    // move the camera forward (in the viewing direction)
    // by the amount dist
	p = p + dist*unit(L);
	setTranslation();
}

void Camera::moveSideways(double dist) {
    // move the camera sideways (orthogonal to the viewing direction)
    // by the amount dist
	p = p + dist*u;
	setTranslation();
}

void Camera::moveVertical(double dist) {
    // move the camera vertically (along the up vector)
    // by the amount dist
	p = p + dist*v;
	setTranslation();
}

void Camera::rotateYaw(double angle) {
    // rotate the camera left/right (around the up vector)
	L = cos(angle) * unit(L) + (-u) * sin(angle);
	setLook(L);
}

void Camera::rotatePitch(double angle) {
    // rotate the camera up/down (pitch angle)
	L = cos(angle) * unit(L) + v*sin(angle);
	v = unit(cross(u, L));
	setRotation();
}

void Camera::rotateAroundAtPoint(int axis, double angle, double focusDist) {
    // Rotate the camera around the right (0), up (1), or look (2) vector
    // around the point at eye + look * focusDist

}


void Camera::moveKeyboard( )
{
    // you may change key controls for the interactive
    // camera controls here, make sure you document your changes
    // in your README file

    if (Fl::event_key('w'))
        moveForward(+0.05);
    if (Fl::event_key('s'))
        moveForward(-0.05);
    if (Fl::event_key('a'))
        moveSideways(-0.05);
    if (Fl::event_key('d'))
        moveSideways(+0.05);
    if (Fl::event_key(FL_Up))
        moveVertical(+0.05);
    if (Fl::event_key(FL_Down))
        moveVertical(-0.05);
    if (Fl::event_key(FL_Left))
        rotateYaw(+0.05);
    if (Fl::event_key(FL_Right))
        rotateYaw(-0.05);
    if (Fl::event_key(FL_Page_Up))
        rotatePitch(+0.05);
    if (Fl::event_key(FL_Page_Down))
        rotatePitch(-0.05);
}

void Camera::setTranslation() {
	T = Matrix4::translation(Point3(-p[0], -p[1], -p[2]));
	T_1 = Matrix4::translation(p);
}

void Camera::setRotation() {
	n = unit(-L);
	v = unit(U - dot(U, n)*n);
	u = unit(cross(v,n));
	R = Matrix4(u, v, n);
	R_1 = R.transpose();
}

void Camera::setFarScaling() {
	Sxyz = Matrix4::scaling(1.0 / df);
	Sxyz_1 = Matrix4::scaling(df);
}

void Camera::setWidthScaling() {
	double ratio = width / 1.0 /height;
	double theta_h = zoom/360.0 * (2 * M_PI);
	Sxy = Matrix4::scaling(1 / (ratio*tan(theta_h / 2.0)), 1 / tan(theta_h / 2.0), 1);
	Sxy_1 = Matrix4::scaling(ratio*tan(theta_h / 2.0), tan(theta_h/2.0), 1);
}

void Camera::setPerspective() {
	double ratio= dn / 1.0 / df;
	D = Matrix4(Vector4(1, 0, 0, 0), Vector4(0, 1, 0, 0), Vector4(0, 0, 1 / (ratio - 1), ratio / (ratio - 1)), Vector4(0, 0, -1, 0));
}