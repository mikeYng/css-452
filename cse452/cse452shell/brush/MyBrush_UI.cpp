#include "cse452.h"
#include "ScreenPoint.h"
#include "BrushInterface.h"
#include <FL/Fl.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/fl_draw.H>
#include <FL/gl.h>
#include <cstring>
#include <cmath>

MyBrush::MyBrush() 
{
    isMouseDown = false;

    imageWidth  = screenWidth = 0;
    imageHeight = screenHeight = 0;

    // initialize your data here
}

MyBrush::~MyBrush() {
    // destroy your data here
}

void MyBrush::resize(int width, int height) {
    screenWidth  = width;
    screenHeight = height;

    // First time initialization
    if ( imageWidth == 0 ) {
        imageWidth = screenWidth;
        imageHeight = screenHeight;

        // Make image black
        pixelData.resize( width * height * 3, 0 );
    }
}

void MyBrush::loadImage(Fl_Image* image) {
    imageWidth = image->w();
    imageHeight = image->h();
    // Reset viewport
    resize( screenWidth, screenHeight );
    pixelData.resize( imageWidth * imageHeight * 3, 0 );

    // OpenGL's windows are reversed in y
    const int delta = imageWidth * 3;
    unsigned char* src = (unsigned char*) *image->data();
    for (int i = 0; i < imageHeight; i++) {
        // Ok, this is ugly
        unsigned char* dest = &pixelData[ ((imageHeight - 1 - i) * imageWidth * 3) ];
        memcpy(dest, src, delta);
        src += delta;
    }
}

void MyBrush::draw() {
    // Set up camera for drawing
    setup2DDrawing( Color(0,0,0), screenWidth, screenHeight );

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Draw a border around the actual image
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i( 0,            0 );
    glVertex2i( imageWidth+1, 0 );
    glVertex2i( imageWidth+1, imageHeight+1 );
    glVertex2i( 0,            imageHeight+1 );
    glEnd();


    glRasterPos2i(0, 0);
    // Copy data into window
	//for ( int iX = 0; iX < 100; iX++ )
		//putPixel( iX, iX, Color(1,0,0) );

    glDrawPixels(imageWidth, imageHeight, GL_RGB, GL_UNSIGNED_BYTE, &pixelData[0]);

	// These 5 lines draw a white line across your canvas
	// Remove this and replace it with intelligent OpenGL preview code

    // Add in your OpenGL pre-view code here

    // display draw in progress (mouse is down)
    if (isMouseDown) {
		int toolType = brushUI->getToolType();
		int radius = brushUI->getRadius();
		Color col = brushUI->getColor();
		switch (toolType) {
		case TOOL_LINE:
			//the UI for the line
			glLineWidth(1.5*radius);
			glColor3f(col[0], col[1], col[2]);
			glBegin(GL_LINES);
			glVertex2i(mouseDown[0], mouseDown[1]);
			glVertex2i(mouseDrag[0], mouseDrag[1]);
			glEnd();
			break;
		case TOOL_POLYGON:
			// The UI for the polygon
			glLineWidth(1);
			glColor3f(col[0], col[1], col[2]);
			glBegin(GL_LINES);
			if (polygon.size() > 0) {
				if (polygon.size() > 1) {
					for (unsigned i = 0; i < polygon.size() - 1; ++i) {
						glVertex2i(polygon[i][0], polygon[i][1]);
						glVertex2i(polygon[i + 1][0], polygon[i + 1][1]);
						//connecting the lines
					}
					//connecting the cursor and the existing points. 
					glVertex2i(polygon[polygon.size() - 1][0], polygon[polygon.size() - 1][1]);
					glVertex2i(mouseDrag[0], mouseDrag[1]);
					glVertex2i(polygon[0][0], polygon[0][1]);
					glVertex2i(mouseDrag[0], mouseDrag[1]);
				}
				else {
					// when only when vertex has been put
					glVertex2i(polygon[0][0], polygon[0][1]);
					glVertex2i(mouseDrag[0], mouseDrag[1]);
				}
			}
			glEnd();
			break;
		case TOOL_CIRCLE:
			// The UI for the circle
			glColor3f(col[0], col[1], col[2]);
			glBegin(GL_LINE_LOOP);
			double d = sqrt(pow(mouseDown[0] - mouseDrag[0], 2) + pow(mouseDown[1] - mouseDrag[1], 2));

			for (int j = (int)d; j < (int)d + radius; j++) {
				// written using polygon approximation. 
				for (int i = 0; i < 100; i++) {
					float angle = 2.0 * 3.1415926 * i / 100.0;
					float x = j * cosf(angle);
					float y = j * sinf(angle);
					glVertex2i(x + mouseDown[0], y + mouseDown[1]);
				}
			}
			glEnd();
			break;
		}
		
    }
	else {
		int toolType = brushUI->getToolType();
		int radius = brushUI->getRadius();
		Color col = brushUI->getColor();
		switch (toolType) {
		case TOOL_BRUSH:
			//UI for the brush. 
			glBegin(GL_LINE_LOOP);
			for (int i = 0; i < 100; i++)
			{
				float angle = 2.0 * 3.1415926 * i / 100.0;

				float x = radius * cosf(angle);
				float y = radius * sinf(angle);
				glVertex2f(x + mouseDrag[0], y + mouseDrag[1]);
			}
			glEnd();
			break;
		}
	}
    endDrawing();
}

// This does pixel flow
void MyBrush::draw_callback( void *in_data )
{
    MyBrush *opMe = static_cast<MyBrush *>( in_data );

    // Repeat the time out if we're not done yet
    if ( opMe->isMouseDown == true ) {
        opMe->drawBrush();

        Fl::repeat_timeout( 0.05, MyBrush::draw_callback, (void *) opMe );

        RedrawWindow();
    }
}


int MyBrush::handle(int event) {
    // OpenGL & FLTK's y axes are oriented differently
    const ScreenPoint pt = ScreenPoint( Fl::event_x(), screenHeight - 1 - Fl::event_y() );

    switch (event) {
        case FL_PUSH: {
            mouseDrag = pt;
            mouseDown = pt;

            if (brushUI->getToolType() == TOOL_POLYGON) {
                if (isMouseDown == true) {
                    polygon.push_back( mouseDrag );
                } else {
                    isMouseDown = true;
                    polygon.resize(0);
                    polygon.push_back( mouseDrag );
                }
            } else {
                isMouseDown = true;
                if (brushUI->getToolType() == TOOL_BRUSH)
                    Fl::add_timeout(0, draw_callback, this);
            }
            return 1;
        }
        case FL_DRAG: mouseDrag = pt; RedrawWindow(); return 1;
        case FL_MOVE: 
            mouseDrag = pt;
            if ( brushUI->getToolType() == TOOL_BRUSH || ( brushUI->getToolType() == TOOL_POLYGON && isMouseDown ) )
                RedrawWindow();
            return 1;
        case FL_RELEASE: {
            mouseDrag = pt;
             if (brushUI->getToolType() != TOOL_POLYGON) {
                isMouseDown = false;
                switch (brushUI->getToolType()) {
                    case TOOL_BRUSH: 
                        break;
                    case TOOL_LINE: 
                        drawLine( ); 
                        break;
                    case TOOL_CIRCLE: 
                        drawCircle( );
                        break;
                    case TOOL_FILTER: 
                        filterRegion( ); 
                        break;
                    default: break;
                }
             } else if ( Fl::event_button3() || Fl::event_state( FL_SHIFT ) ) {
                 isMouseDown = false;
                 drawPolygon();
             }
             RedrawWindow();
            return 1;
        }
        default: return 0;
    }
}
