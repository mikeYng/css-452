#include "cse452.h"
#include "MyBrush.h"
#include "BrushInterface.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

void MyBrush::PutPixel(int x, int y, Color col) {
	//the wrapper for putting pixel, with boundary checks
	if (x > 0 && x < screenWidth && y > 0 && y < screenHeight) {
		putPixel(x, y, col);
	}
}

int MyBrush::maxPointY(vector<ScreenPoint> a) {
	//Getting the max Y value from vector of points. 
	ScreenPoint maxY;
	maxY[1] = 0;
	for (unsigned i = 0; i < a.size(); ++i) {
		if (a[i][1] > maxY[1]) {
			maxY = a[i];
		}
	}
	return maxY[1];
}

int MyBrush::maxPointX(vector<ScreenPoint> a) {
	//Getting the max X value from vector of points
	ScreenPoint maxX;
	maxX[1] = 0;
	for (unsigned i = 0; i < a.size(); ++i) {
		if (a[i][0] > maxX[0]) {
			maxX = a[i];
		}
	}
	return maxX[0];
}

int MyBrush::minPointX(vector<ScreenPoint> a) {
	//Getting the min X value from vector of points
	ScreenPoint minX;
	minX[0] = screenWidth + 1;
	for (unsigned i = 0; i < a.size(); ++i) {
		if (a[i][0] < minX[0]) {
			minX = a[i];
		}
	}
	return minX[0];
}

int MyBrush::minPointY(vector<ScreenPoint> a) {
	//Getting the min Y value from vector of points
	ScreenPoint minY;
	minY[1] = screenHeight+1;
	for (unsigned i = 0; i < a.size(); ++i) {
		if (a[i][1] < minY[1]) {
			minY = a[i];
		}
	}
	return minY[1];
}
int MyBrush::scanLine(int y, double slope, ScreenPoint point) {
	//getting the intersections with some line. 
	return (int)round((y - point[1]) / slope + point[0]);

}

void MyBrush::horizontal(int radius, Color col) {
	//horizontal case
	int upperY = mouseDown[1] + radius;
	int lowerY = mouseDown[1] - radius;
	for (int i = 0; i < 2 * radius; ++i) {
		for (int j = 0; j < abs(mouseDown[0] - mouseDrag[0]); ++j) {
			if (mouseDrag[0] > mouseDown[0]) {
				PutPixel(mouseDown[0] + j, lowerY + i, col);
			}
			else {
				PutPixel(mouseDown[0] - j, lowerY + i, col);
			}

		}
	}
}

void MyBrush::vertical(int radius, Color col) {
	//vertical case for the line. 
	int leftX = mouseDown[0] - radius;
	int rightX = mouseDown[0] + radius;
	for (int i = 0; i < 2 * radius; ++i) {
		for (int j = 0; j < abs(mouseDown[1] - mouseDrag[1]); ++j) {
			if (mouseDrag[1] > mouseDown[1]) {
				PutPixel(leftX + i, mouseDown[1] + j, col);
			}
			else {
				PutPixel(leftX + i, mouseDown[1] - j, col);
			}

		}
	}
}
int MyBrush::isvertex(int y, vector<ScreenPoint>points, int j) {
	for (unsigned i = 0; i < points.size(); ++i) {
		if (points[i][1] == y && (i==j || i ==j+1)) {
			return points[i][0];
		}
	}
	return -1;
}
void MyBrush::scanRegion(int maxY, int minY, vector<ScreenPoint>points, Color col) {
	// function i use to fill in polygons. 
	vector<int> intersections;
	for (int i = 0; i <= maxY - minY; ++i) {
		for (unsigned j = 0; j < points.size() - 1; ++j) {
			if ((minY + i >= points[j][1] && minY + i < points[j + 1][1]) || (minY + i < points[j][1] && minY + i >= points[j + 1][1]) || (points[j][0] == points[j + 1][0])) {
				//if the scan line is in between that edge. 
				if (points[j][1] == points[j + 1][1]) {
					intersections.push_back(points[j][0]);
					intersections.push_back(points[j + 1][0]);
				}
				else if (points[j][0] == points[j + 1][0]) {
					intersections.push_back(points[j][0]);
				}
				else {
					intersections.push_back(scanLine(minY + i, (points[j + 1][1] - points[j][1]) / 1.0 / (points[j + 1][0] - points[j][0]), points[j]));
				}

			}
		}
		sort(intersections.begin(), intersections.end());
		//intersections.erase(unique(intersections.begin(), intersections.end()), intersections.end());
		if (intersections.size() >= 1) {
			for (unsigned j = 0; j < intersections.size() - 1; j += 2) {
				//since we only want to draw things from odd intersections to even ones. 
				if (j + 1 < intersections.size()) {
					int spread = abs(intersections[j] - intersections[j + 1]);
					int left;
					if (intersections[j] < intersections[j + 1]) {
						left = intersections[j];
					}
					else {
						left = intersections[j + 1];
					}
					for (int k = 0; k <= spread; ++k) {
						PutPixel(left + k, minY + i, col);
					}
				}

			}
			intersections.clear();
		}

	}
}
void MyBrush::changedBrush() {
	// this is called anytime the brush type or brush radius changes
	// it should recompute the brush mask appropriately
	const int radius = brushUI->getRadius();
	const int brushType = brushUI->getBrushType();
	mask = vector<vector<double>>(); //need (radius+1)^2 mask pixel data
	mask.resize(radius + 1, vector<double>(radius + 1, 0));
	switch (brushType) {
	case BRUSH_CONSTANT:
		for (int i = 0; i <= radius; ++i) {
			for (int j = 0; j <= radius; ++j) {
				if (i*i + j*j <= radius*radius) {
					//within the circle
					mask[i][j] = 1;
				}
				else {
					//outside the circle
					mask[i][j] = 0;
				}
			}
		}
		break;
	case BRUSH_LINEAR:
		for (int i = 0; i <= radius; ++i) {
			for (int j = 0; j <= radius; ++j) {
				if (i*i + j*j <= radius*radius) {
					double d = sqrt(i*i + j*j);
					mask[i][j] = (radius - d) /1.0/ radius;
				}
				else {
					mask[i][j] = 0;
				}
			}
		}
		break;
	case BRUSH_QUADRATIC:
		for (int i = 0; i <= radius; ++i) {
			for (int j = 0; j <= radius; ++j) {
				if (i*i + j*j <= radius*radius) {
					double d = sqrt(i*i + j*j);
					double ratio = (d*1.0) / radius;
					mask[i][j] = pow(1-ratio,2);
				}
				else {
					mask[i][j] = 0;
				}
			}
		}
		break;
	case BRUSH_GAUSSIAN:
		for (int i = 0; i <= radius; ++i) {
			for (int j = 0; j <= radius; ++j) {
				if (i*i + j*j <= radius*radius) {
					//within the circle
					double d = sqrt(i*i + j*j);
					mask[i][j] = 2.5/sqrt(2*M_PI)*exp(-(pow(2.5*d/radius,2))/2);
				}
				else {
					//outside the circle
					mask[i][j] = 0;
				}
			}
		}
		break;
	case BRUSH_SPECIAL:
		for (int i = 0; i <= radius; ++i) {
			for (int j = 0; j <= radius; ++j) {
				if (i*i + j*j <= radius*radius) {
					//within the circle
					double d = sqrt(i*i + j*j);
					double ratio = (d*1.0) / radius;
					mask[i][j] = 1.0-pow(1 - ratio, 2);
				}
				else {
					//outside the circle
					mask[i][j] = 0;
				}
			}
		}
	}
}

void MyBrush::drawBrush( ) {
    // apply the current brush mask to image location (x,y)
    // the mouse location is in mouseDrag

    const int radius = brushUI->getRadius();
    const float pixelFlow = brushUI->getFlow();
    const Color colBrush = brushUI->getColor();
	Color actual;
	Color canvas;
	for (int i = 0; i <=radius; ++i) {
		for (int j = 0; j <= radius; ++j) {
			if (mask[i][j] != 0) {
				if (i != 0 && j != 0) {
					// normal points and checking the boundaries. 
					if ((mouseDrag[0] + i < screenWidth) && (mouseDrag[1] + j < screenHeight)) {
						canvas = getPixel(mouseDrag[0] + i, mouseDrag[1] + j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;//pixel flow 
						putPixel(mouseDrag[0] + i, mouseDrag[1] + j, actual);

					}
					if ((mouseDrag[0] - i >= 0) && (mouseDrag[1] + j < screenHeight)) {
						canvas = getPixel(mouseDrag[0] - i, mouseDrag[1] + j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] - i, mouseDrag[1] + j, actual);
					}
					if ((mouseDrag[0] + i < screenWidth) && (mouseDrag[1] - j >= 0)) {
						canvas = getPixel(mouseDrag[0] + i, mouseDrag[1] - j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] + i, mouseDrag[1] - j, actual);
					}
					if ((mouseDrag[0] - i >= 0) && (mouseDrag[1] - j >= 0)) {
						canvas = getPixel(mouseDrag[0] - i, mouseDrag[1] - j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] - i, mouseDrag[1] - j, actual);
					}
				}
				else if (i == 0 && j != 0) {
					//these are the two cross lines. 
					if (mouseDrag[1] + j < screenHeight) {
						canvas = getPixel(mouseDrag[0] , mouseDrag[1] + j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] , mouseDrag[1] + j, actual);
					}
					if (mouseDrag[1] - j >= 0) {
						canvas = getPixel(mouseDrag[0], mouseDrag[1] - j);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0], mouseDrag[1] - j, actual);
					}
				}
				else if (i != 0 && j == 0) {
					if (mouseDrag[0] - i >= 0) {
						canvas = getPixel(mouseDrag[0] - i, mouseDrag[1]);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] - i, mouseDrag[1], actual);
					}
					if (mouseDrag[0] + i < screenWidth) {
						canvas = getPixel(mouseDrag[0] + i, mouseDrag[1]);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0] + i, mouseDrag[1], actual);
					}
				}
				else {
					// the center 
					if (mouseDrag[1] < screenHeight && mouseDrag[0] < screenWidth && mouseDrag[0] > 0 && mouseDrag[1] > 0 ) {
						canvas = getPixel(mouseDrag[0], mouseDrag[1]);
						actual = colBrush*mask[i][j] * pixelFlow + (1 - mask[i][j] * pixelFlow)*canvas;
						putPixel(mouseDrag[0], mouseDrag[1], actual);
					}
					
				}
				
			}
		}
	}
}

void MyBrush::drawLine( ) {
    // draw a thick line from mouseDown to mouseDrag
    // the width of the line is given by the current brush radius
    const int radius = brushUI->getRadius();
    const Color colBrush = brushUI->getColor();
	if (mouseDrag[0] != mouseDown[0]) {
		double slope = (mouseDrag[1] - mouseDown[1]) /1.0/ (mouseDrag[0] - mouseDown[0]); //getting the slope btw the the two points
		//cout << mouseDrag[0] << ", " << mouseDrag[1] << " " << mouseDown[0] <<", " << mouseDown[1]<< endl;
		if (slope == 0) {
			horizontal(radius, colBrush);
		}
		else {
			double perpenSlope = -1.0 / slope; //get the slope of line perpendicular to it.
			double angle = atan(perpenSlope);
			//cout << "slope: " << slope << " perpen: " << perpenSlope << endl;
			ScreenPoint point1;
			ScreenPoint point2;
			ScreenPoint point3;
			ScreenPoint point4;
			point1[0] = mouseDown[0] + round(cos(angle)*radius / 2.0);
			point1[1] = mouseDown[1] + round(sin(angle)*radius / 2.0);
			point2[0] = mouseDown[0] - round(cos(angle)*radius / 2.0);
			point2[1] = mouseDown[1] - round(sin(angle)*radius / 2.0);
			point3[0] = mouseDrag[0] - round(cos(angle)*radius / 2.0);
			point3[1] = mouseDrag[1] - round(sin(angle)*radius / 2.0);
			point4[0] = mouseDrag[0] + round(cos(angle)*radius / 2.0);
			point4[1] = mouseDrag[1] + round(sin(angle)*radius / 2.0);
			//the commented lines below are equivalent but somehow lack a little accuracy in c++. 

			//point1[0] = (int)round(mouseDown[0]*1.0 - radius / sqrt(1 + pow(perpenSlope, 2))); //get points that are on that line but r distance away from mouseDrag. 
			//point1[1] = (int)round(perpenSlope*(point1[0] - mouseDown[0]*1.0) + mouseDown[1]);
			//point2[0] = (int)round(mouseDown[0]*1.0 + radius / sqrt(1 + pow(perpenSlope, 2)));
			//point2[1] = (int)round(perpenSlope*(point2[0] - mouseDown[0]*1.0) + mouseDown[1]);
			//point3[0] = (int)round(mouseDrag[0]*1.0 + radius / sqrt(1 + pow(perpenSlope, 2)));
			//point3[1] = (int)round(perpenSlope*(point3[0] - mouseDrag[0]*1.0) + mouseDrag[1]);
			//point4[0] = (int)round(mouseDrag[0]*1.0 - radius / sqrt(1 + pow(perpenSlope, 2)));
			//point4[1] = (int)round(perpenSlope*(point4[0] - mouseDrag[0]*1.0) + mouseDrag[1]);
			//putting points in the order such that point i and point i+1 will always have an edge btw them. 
			vector<ScreenPoint> points;
			points.push_back(point1);
			points.push_back(point2);
			points.push_back(point3);
			points.push_back(point4);
			points.push_back(point1);
			int maxY = maxPointY(points);
			int minY = minPointY(points);
			vector<int> intersections;
			for (int i = 0; i <= maxY - minY; ++i) {
				for (unsigned j = 0; j < points.size() - 1; ++j) {
					if ((minY + i >= points[j][1] && minY + i <= points[j + 1][1]) || (minY + i <= points[j][1] && minY + i >= points[j + 1][1])) {
						//if the scan line is in between that edge. 
						if (points[j][1] == points[j + 1][1]) {
							intersections.push_back(points[j][0]);
							intersections.push_back(points[j + 1][0]);
						}
						else if (points[j][0] == points[j + 1][0]) {
							intersections.push_back(points[j][0]);
						}
						else {
							intersections.push_back(scanLine(minY + i, (points[j + 1][1] - points[j][1]) / 1.0 / (points[j + 1][0] - points[j][0]), points[j]));
						}

					}
				}
				sort(intersections.begin(), intersections.end());
				intersections.erase(unique(intersections.begin(), intersections.end()), intersections.end());
				//cout << "size before: " << intersections.size() << endl;
				if (intersections.size() >= 1) {
					for (unsigned j = 0; j < intersections.size() - 1; j += 2) {
						//since we only want to draw things from odd intersections to even ones. 
						if (j + 1 < intersections.size()) {
							int spread = abs(intersections[j] - intersections[j + 1]);
							int left;
							if (intersections[j] < intersections[j + 1]) {
								left = intersections[j];
							}
							else {
								left = intersections[j + 1];
							}
							for (int k = 0; k <= spread; ++k) {
								if (left + k >= minPointX(points) && left + k <= maxPointX(points) && minY + i <= maxY) {
									//cout << "drawing?!" << endl;
									PutPixel(left + k, minY + i, colBrush);
								}
							}
						}

					}
					intersections.clear();
				}

			}

		}
	}
	else {
		vertical(radius, colBrush);
	}
	
	
}


void MyBrush::drawCircle() {
    // draw a thick circle at mouseDown with radius r
    // the width of the circle is given by the current brush radius
	const int radius = brushUI->getRadius();
	const Color colBrush = brushUI->getColor();
	double dis = sqrt(pow(mouseDown[0] - mouseDrag[0], 2) + pow(mouseDown[1] - mouseDrag[1], 2));
	int minY = (int)round(mouseDown[1] - radius - dis);
	int maxY = (int)round(mouseDown[1] + radius + dis);
	int innerDown = (int)round(mouseDown[1] - dis);
	int innerUp = (int)round(mouseDown[1] + dis);
	vector<int> intersections;
	for (int i = minY; i <= maxY ; ++i) {
		// case for 2 intersections. 
		if ((i > minY && i < innerDown) || (i >innerUp && i<maxY)) {
			intersections.push_back((int)round(-sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back((int)round(sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
		}
		else if (i == minY || i == maxY) {
			intersections.push_back(mouseDown[0]);
			intersections.push_back(mouseDown[0]);
		}
		else if (i == innerDown || i==innerUp) {
			// case for 3 intersections. 
			intersections.push_back((int)round(-sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back((int)round(sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back(mouseDown[0]);
			intersections.push_back(mouseDown[0]);
		}
		else {
			//case for 4 intersections. 
			intersections.push_back((int)round(-sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back((int)round(sqrt(pow(dis + radius, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back((int)round(-sqrt(pow(dis, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
			intersections.push_back((int)round(sqrt(pow(dis, 2) - pow(i - mouseDown[1], 2)) + mouseDown[0]));
		}
		sort(intersections.begin(), intersections.end());
		if (intersections.size() >= 1) {
			for (unsigned j = 0; j < intersections.size() - 1; j += 2) {
				//since we only want to draw things from odd intersections to even ones. 
				if (j + 1 < intersections.size()) {
					int spread = abs(intersections[j] - intersections[j + 1]);
					int left;
					if (intersections[j] < intersections[j + 1]) {
						left = intersections[j];
					}
					else {
						left = intersections[j + 1];
					}
					for (int k = 0; k <= spread; ++k) {
						PutPixel(left + k, i, colBrush);
					}
				}

			}
			intersections.clear();
		}

	}
}


void MyBrush::drawPolygon() {
    // draw a polygon with numVertices whos coordinates are stored in the
    // polygon array: {x0, y0, x1, y1, ...., xn-1, yn-1}
	const Color colBrush = brushUI->getColor();
	polygon.push_back(polygon[0]);
	int maxY = maxPointY(polygon);
	int minY = minPointY(polygon);
	scanRegion(maxY, minY, polygon, colBrush);
}

void MyBrush::filterRegion( ) {
    // apply the filter indicated by filterType to the square
    // defined by the two corner points mouseDown and mouseDrag
    // these corners are not guarenteed to be in any order
    // The filter width is given by the brush radius
}
