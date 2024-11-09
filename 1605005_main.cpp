#include <windows.h>
#include <GL/glut.h>

#include <bits/stdc++.h>

#include "1605005_classes.h"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

#define L 0
#define U 1
#define R 2
#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3
#define COUNTER_CLOCKWISE 4
#define CLOCKWISE 5

using namespace std;

point l = {sqrt(.5), sqrt(0.5), 0}, u = {0, 0, 1}, r = {sqrt(.5), -sqrt(0.5), 0}, camera = {-100, -100, 0};

double cameraHeight;
double cameraAngle;
bool drawaxes = true;
double k = 3, rotation_angle = 3.0*pi/180.0;

int recursion_level, pixels, object_count;

vector<Object*> objects;
vector<Light*> lights;

double toRadian(double deg) {
    return deg*pi/180.0;
}

void drawAxes()
{
	if (drawaxes) {
		glColor3f(1.0, 1.0, 1.0);

		glBegin(GL_LINES); {
			glVertex3f( 1000,0,0);
			glVertex3f(-1000,0,0);

			glVertex3f(0,-1000,0);
			glVertex3f(0, 1000,0);

			glVertex3f(0,0, 3000);
			glVertex3f(0,0,-3000);

		} glEnd();
	}
}

void rotateCamera(int axis, double angle) {
    // rotate the other two axes w.r.t the specified axis

    point t1, t2;
    if (axis == U) {
        // look left/right
        // rotate L first
        t1 = l*cos(angle) + r*sin(angle);
        // rotate R now
        t2 = r*cos(angle) - l*sin(angle);

        l = t1;
        r = t2;
    }
    else if (axis == R) {
        // look up/down
        // rotate L first
        t1 = l*cos(angle) + u*sin(angle);
        // rotate U now
        t2 = u*cos(angle) - l*sin(angle);

        l = t1;
        u = t2;
    }
    else if (axis == L) {
        // tilt clockwise/counter-clockwise
       // rotate R first
        t1 = r*cos(angle) + u*sin(angle);
        // rotate U now
        t2 = u*cos(angle) - r*sin(angle);

        r = t1;
        u = t2;
    }
}

void look(int direction) {
    if (direction == LEFT) {
        rotateCamera(U, -rotation_angle);
    }
    else if(direction == RIGHT) {
        rotateCamera(U, rotation_angle);
    }
    else if (direction == UP) {
        rotateCamera(R, rotation_angle);
    }
    else if (direction == DOWN) {
        rotateCamera(R, -rotation_angle);
    }
    else if (direction == COUNTER_CLOCKWISE) {
        rotateCamera(L, -rotation_angle);
    }
    else if (direction == CLOCKWISE) {
        rotateCamera(L, rotation_angle);
    }
}

void draw() {
    for (int i=0; i<objects.size(); i++) {
        objects[i]->draw();
	}

	for (int i=0; i<lights.size(); i++) {
        lights[i]->draw();
	}
}

void capture() {
    cout << "Capturing" << endl;
    double imageWidth = pixels, imageHeight = pixels;
    bitmap_image image(imageWidth, imageHeight);
    for (int i=0; i<imageWidth; i++) {
        for (int j=0; j<imageHeight; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    point eye = camera;
    double windowHeight = 500, windowWidth = 500, viewAngle = toRadian(80);
    double planeDistance = (windowHeight/2.0) * tan(viewAngle/2.0);
    point topLeft = eye + l*planeDistance - r*(windowWidth/2.0) + u*(windowHeight/2.0);

    double du = windowWidth/imageWidth;
    double dv = windowHeight/imageHeight;
    topLeft = topLeft + r*(du/2.0) - u*(dv/2.0);

    int nearest;
    double t, tMin = 10000000;
    Object* oMin = NULL;
    Object* o;
    for (int i=0; i<imageWidth; i++) {
        for (int j=0; j<imageHeight; j++) {
            point curPixel = topLeft + r*(i*du) - u*(j*dv);
            Ray ray(eye, curPixel-eye);
            double* color = new double[3];
            tMin = 1000000;
            oMin = 0;
            for (int k=0; k<objects.size(); k++) {
                o = objects[k];
                t = o->intersect(&ray, color, 0);
                if (t > 0) {
                    if (t < tMin) {
                        tMin = t;
                        oMin = objects[k];
                    }
                }
            }

            if (oMin) {
                t = oMin->intersect(&ray, color, 1);
                color[0] = color[0] < 0 ? 0 : color[0] > 1 ? 1 : color[0];
                color[1] = color[1] < 0 ? 0 : color[1] > 1 ? 1 : color[1];
                color[2] = color[2] < 0 ? 0 : color[2] > 1 ? 1 : color[2];
                image.set_pixel(i, j, 255*color[0], 255*color[1], 255*color[2]);
            }

            delete[] color;
            color = 0;
        }
    }

//    image.save_image("H:\\Academic\\4-1\\output.bmp");
    image.save_image("H:\\Academic\\4-1\\CSE 409 (Graphics)\\Lab\\Offline 3\\Ray Tracing\\Ray Tracing\\1605005_output.bmp");
    cout << "Saved" << endl;
}

void keyboardListener(unsigned char key, int x, int y) {
    switch(key) {
    case '1':
        look(LEFT);
        break;
    case '2':
        look(RIGHT);
        break;
    case '3':
        look(UP);
        break;
    case '4':
        look(DOWN);
        break;
    case '5':
        look(CLOCKWISE);
        break;
    case '6':
        look(COUNTER_CLOCKWISE);
        break;
    case '0':
        capture();
        break;
    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y) {
	switch(key) {
		case GLUT_KEY_UP:		// up arrow key
			// move forward
			camera.x += k * l.x;
			camera.y += k * l.y;
			camera.z += k * l.z;
			break;

        case GLUT_KEY_DOWN:		//down arrow key
            // move backward
			camera.x -= k * l.x;
			camera.y -= k * l.y;
			camera.z -= k * l.z;
			break;

		case GLUT_KEY_RIGHT:
			// move right
			camera.x += k * r.x;
			camera.y += k * r.y;
			camera.z += k * r.z;
			break;

		case GLUT_KEY_LEFT:
			// move left
			camera.x -= k * r.x;
			camera.y -= k * r.y;
			camera.z -= k * r.z;
			break;

		case GLUT_KEY_PAGE_UP:
			// move up
			camera.x += k * u.x;
			camera.y += k * u.y;
			camera.z += k * u.z;
			break;

		case GLUT_KEY_PAGE_DOWN:
			// move left
			camera.x -= k * u.x;
			camera.y -= k * u.y;
			camera.z -= k * u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y) {
	switch(button) {
		case GLUT_LEFT_BUTTON:
		    if (state == GLUT_DOWN) {

		    }
			break;

		case GLUT_RIGHT_BUTTON:
		    if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes ^= 1;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display() {

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?
	gluLookAt(camera.x, camera.y, camera.z,	camera.x+l.x, camera.y+l.y, camera.z+l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	draw();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate() {
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init() {
	//codes for initialization
	drawaxes = 1;
	cameraHeight = 150.0;
	cameraAngle = 1.0;

	//clear the screen
	glClearColor(0, 0, 0, 0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80, 1, 1, 1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void loadData() {
    ifstream ifp;
//    ifp.open("H:\\Academic\\4-1\\scene.txt");
    ifp.open("H:\\Academic\\4-1\\CSE 409 (Graphics)\\Lab\\Offline 3\\Ray Tracing\\Ray Tracing\\scene.txt");
    ifp >> recursion_level >> pixels >> object_count;

    string obj;
    Object* temp;
    point p1, p2, p3;
    double v1;
    for (int i=0; i<object_count; i++) {
        ifp >> obj;
        if (obj == "sphere") {
            ifp >> p1.x >> p1.y >> p1.z;
            ifp >> v1;
            temp = new Sphere(p1, v1);
            ifp >> temp->color[0] >> temp->color[1] >> temp->color[2];
            ifp >> temp->coEfficients[0] >> temp->coEfficients[1] >> temp->coEfficients[2] >> temp->coEfficients[3];
            ifp >> temp->shine;
            objects.push_back(temp);
        } else if (obj == "triangle") {
            ifp >> p1.x >> p1.y >> p1.z;
            ifp >> p2.x >> p2.y >> p2.z;
            ifp >> p3.x >> p3.y >> p3.z;
            temp = new Triangle(p1, p2, p3);
            ifp >> temp->color[0] >> temp->color[1] >> temp->color[2];
            ifp >> temp->coEfficients[0] >> temp->coEfficients[1] >> temp->coEfficients[2] >> temp->coEfficients[3];
            ifp >> temp->shine;
            objects.push_back(temp);
        } else if (obj == "general") {
            temp = new General();
            ifp >> temp->a >> temp->b >> temp->c >> temp->d >> temp->e >> temp->f >> temp->g >> temp->h >> temp->i >> temp->j;
            ifp >> p1.x >> p1.y >> p1.z;
            temp->reference_point = p1;
            ifp >> temp->length >> temp->width >> temp->height;
            ifp >> temp->color[0] >> temp->color[1] >> temp->color[2];
            ifp >> temp->coEfficients[0] >> temp->coEfficients[1] >> temp->coEfficients[2] >> temp->coEfficients[3];
            ifp >> temp->shine;
            objects.push_back(temp);
        }
    }
    temp = new Floor(1000, 20);
    temp->coEfficients[0] = 0.4;
    temp->coEfficients[1] = 0.2;
    temp->coEfficients[2] = 0.1;
    temp->coEfficients[3] = 0.2;
    temp->shine = 5;
    objects.push_back(temp);

    int light_count;
    ifp >> light_count;
    for (int i=0; i<light_count; i++) {
        Light* light = new Light();
        ifp >> light->pos.x >> light->pos.y >> light->pos.z;
        ifp >> light->color[0] >> light->color[1] >> light->color[2];
        lights.push_back(light);
    }

    ifp.close();
}


void clearVectors() {
    cout << "Clearing" << endl;
    for (int i=0; i<objects.size(); i++) {
        Object* o = objects[i];
        delete o;
    }

    for (int i=0; i<lights.size(); i++) {
        Light* l = lights[i];
        delete l;
    }

    cout << "Cleared" << endl;
}

int main(int argc, char **argv) {
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracer");

	loadData();
	atexit(clearVectors);
	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occurring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
