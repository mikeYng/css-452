// generated by Fast Light User Interface Designer (fluid) version 1.0107

#include "IntersectionInterface.h"

void IntersectionInterface::cb_Sphere_i(Fl_Menu_*, void*) {
  intersectionUI.changeShape( (ShapesUI::ShapeType) m_iShapeType->value());
RedrawWindow();
}
void IntersectionInterface::cb_Sphere(Fl_Menu_* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_Sphere_i(o,v);
}

void IntersectionInterface::cb_Cone_i(Fl_Menu_*, void*) {
  intersectionUI.changeShape( (ShapesUI::ShapeType) m_iShapeType->value());
RedrawWindow();
}
void IntersectionInterface::cb_Cone(Fl_Menu_* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_Cone_i(o,v);
}

void IntersectionInterface::cb_Cylinder_i(Fl_Menu_*, void*) {
  intersectionUI.changeShape( (ShapesUI::ShapeType) m_iShapeType->value());
RedrawWindow();
}
void IntersectionInterface::cb_Cylinder(Fl_Menu_* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_Cylinder_i(o,v);
}

void IntersectionInterface::cb_Cube_i(Fl_Menu_*, void*) {
  intersectionUI.changeShape( (ShapesUI::ShapeType) m_iShapeType->value());
RedrawWindow();
}
void IntersectionInterface::cb_Cube(Fl_Menu_* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_Cube_i(o,v);
}

Fl_Menu_Item IntersectionInterface::menu_m_iShapeType[] = {
 {"Sphere", 0,  (Fl_Callback*)IntersectionInterface::cb_Sphere, 0, 4, FL_NORMAL_LABEL, 0, 14, 0},
 {"Cone", 0,  (Fl_Callback*)IntersectionInterface::cb_Cone, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {"Cylinder", 0,  (Fl_Callback*)IntersectionInterface::cb_Cylinder, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {"Cube", 0,  (Fl_Callback*)IntersectionInterface::cb_Cube, 0, 0, FL_NORMAL_LABEL, 0, 14, 0},
 {0,0,0,0,0,0,0,0,0}
};

void IntersectionInterface::cb_m_dXAt_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dXAt(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dXAt_i(o,v);
}

void IntersectionInterface::cb_m_dYAt_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dYAt(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dYAt_i(o,v);
}

void IntersectionInterface::cb_m_dZAt_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dZAt(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dZAt_i(o,v);
}

void IntersectionInterface::cb_m_dTheta_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dTheta(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dTheta_i(o,v);
}

void IntersectionInterface::cb_m_dPhi_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dPhi(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dPhi_i(o,v);
}

void IntersectionInterface::cb_Write_i(Fl_Button*, void*) {
  intersectionUI.writeTest();
}
void IntersectionInterface::cb_Write(Fl_Button* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_Write_i(o,v);
}

void IntersectionInterface::cb_m_dXRot_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dXRot(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dXRot_i(o,v);
}

void IntersectionInterface::cb_m_dYRot_i(Fl_Value_Slider*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_dYRot(Fl_Value_Slider* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_dYRot_i(o,v);
}

void IntersectionInterface::cb_m_bGrid_i(Fl_Check_Button*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_bGrid(Fl_Check_Button* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_bGrid_i(o,v);
}

void IntersectionInterface::cb_m_bRay_i(Fl_Check_Button*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_bRay(Fl_Check_Button* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_bRay_i(o,v);
}

void IntersectionInterface::cb_m_bRayShadow_i(Fl_Check_Button*, void*) {
  RedrawWindow();
}
void IntersectionInterface::cb_m_bRayShadow(Fl_Check_Button* o, void* v) {
  ((IntersectionInterface*)(o->parent()->user_data()))->cb_m_bRayShadow_i(o,v);
}

Fl_Double_Window* IntersectionInterface::make_window() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = m_intersectionWindow = new Fl_Double_Window(420, 265, "Intersection UI");
    w = o;
    o->user_data((void*)(this));
    { Fl_Group* o = new Fl_Group(5, 25, 145, 30);
      o->end();
    }
    { Fl_Choice* o = m_iShapeType = new Fl_Choice(5, 25, 145, 30, "Object type");
      o->down_box(FL_BORDER_BOX);
      o->align(FL_ALIGN_TOP_LEFT);
      o->menu(menu_m_iShapeType);
    }
    { Fl_Value_Slider* o = m_dXAt = new Fl_Value_Slider(5, 75, 200, 25, "At x pos");
      o->type(5);
      o->minimum(-1.5);
      o->maximum(1.5);
      o->callback((Fl_Callback*)cb_m_dXAt);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = m_dYAt = new Fl_Value_Slider(5, 115, 200, 25, "At y pos");
      o->type(5);
      o->minimum(-1.5);
      o->maximum(1.5);
      o->callback((Fl_Callback*)cb_m_dYAt);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = m_dZAt = new Fl_Value_Slider(5, 155, 200, 25, "At z pos");
      o->type(5);
      o->minimum(-1.5);
      o->maximum(1.5);
      o->callback((Fl_Callback*)cb_m_dZAt);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = m_dTheta = new Fl_Value_Slider(5, 195, 200, 25, "Vec theta");
      o->type(5);
      o->maximum(360);
      o->step(1);
      o->callback((Fl_Callback*)cb_m_dTheta);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = m_dPhi = new Fl_Value_Slider(5, 235, 200, 25, "Vec phi");
      o->type(5);
      o->minimum(-90);
      o->maximum(90);
      o->step(1);
      o->value(45);
      o->callback((Fl_Callback*)cb_m_dPhi);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Button* o = new Fl_Button(330, 25, 85, 25, "Write test");
      o->callback((Fl_Callback*)cb_Write);
    }
    { Fl_Value_Slider* o = m_dXRot = new Fl_Value_Slider(215, 75, 200, 25, "View rotation");
      o->type(5);
      o->maximum(360);
      o->step(1);
      o->callback((Fl_Callback*)cb_m_dXRot);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Value_Slider* o = m_dYRot = new Fl_Value_Slider(215, 115, 200, 25, "View height");
      o->type(5);
      o->minimum(-90);
      o->maximum(90);
      o->step(1);
      o->callback((Fl_Callback*)cb_m_dYRot);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    { Fl_Check_Button* o = m_bGrid = new Fl_Check_Button(215, 155, 25, 25, "Show grid");
      o->down_box(FL_DOWN_BOX);
      o->value(1);
      o->callback((Fl_Callback*)cb_m_bGrid);
    }
    { Fl_Check_Button* o = m_bRay = new Fl_Check_Button(215, 195, 25, 25, "Show ray");
      o->down_box(FL_DOWN_BOX);
      o->value(1);
      o->callback((Fl_Callback*)cb_m_bRay);
    }
    { Fl_Check_Button* o = m_bRayShadow = new Fl_Check_Button(215, 235, 25, 25, "Show ray shadow");
      o->down_box(FL_DOWN_BOX);
      o->value(1);
      o->callback((Fl_Callback*)cb_m_bRayShadow);
    }
    m_iSeed = new Fl_Value_Input(240, 30, 85, 20, "Seed");
    o->end();
    o->resizable(o);
  }
  return w;
}

IntersectionInterface::IntersectionInterface() {
  m_intersectionWindow = make_window();
intersectionUI.setUI( this );
intersectionUI.changeShape( ShapesUI::SHAPE_SPHERE );
}

UIInterface * IntersectionInterface::getUI() {
  return &intersectionUI;
}

double IntersectionInterface::getTheta() const {
  return M_PI * m_dTheta->value() / 180.0;
}

double IntersectionInterface::getPhi() const {
  return M_PI * m_dPhi->value() / 180.0;
}

double IntersectionInterface::getXRot() const {
  return M_PI * m_dXRot->value() / 180.0;
}

double IntersectionInterface::getYRot() const {
  return M_PI * m_dYRot->value() / 180.0;
}
