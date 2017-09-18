#ifndef _INTERSECTION_UI_H_
#define _INTERSECTION_UI_H_

#include "../cse452.h"
#include "../shapes/ShapesUI.h"
#include "HitRecord.h"
#include "../UIInterface.h"
#include <FL/Fl_Window.H>
#include <memory>
using namespace std;
class IntersectionInterface;
class IntersectionUI : public UIInterface {
public:
    IntersectionUI();
    ~IntersectionUI();

    // Inherited from userInterface
    void resize(int width, int height);
    void draw();
    int handle(int event);

    // Link to the intersection widget
    void setUI( const IntersectionInterface *in_ui ) { intersectionUI = in_ui; }
    void changeShape( ShapesUI::ShapeType type );
	HitRecord intersect(Point3, Vector3);
    void writeTest();

private:
    const IntersectionInterface *intersectionUI;
    int width, height;

    void drawHits(HitRecord& hr);
	shared_ptr<Shape> shape;
    // declare your variables here
};

#endif /* _INTERSECTION_UI_H_ */
