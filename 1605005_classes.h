#include <bits/stdc++.h>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
const double eps = 1e-4;
extern int recursion_level;

using namespace std;

struct point
{
	double x, y, z;

	point() {
	    x = y = z = 0;
	}

	point(double x, double y, double z) {
	    this->x = x;
	    this->y = y;
	    this->z = z;
	}

	point operator+(point b) {
	    return {x+b.x, y+b.y, z+b.z};
	}

    point operator-(point b) {
        return {x-b.x, y-b.y, z-b.z};
    }

	point operator-() {
	    return {-x, -y, -z};
	}

	point operator*(double scalar) {
	    return {x*scalar, y*scalar, z*scalar};
	}

	point operator*(point& rhs) {
        // cross
        return point(y*rhs.z-z*rhs.y, z*rhs.x-x*rhs.z, x*rhs.y-y*rhs.x);
    }

    bool operator==(point& rhs) {
        return abs(x-rhs.x)<eps && abs(y-rhs.y)<eps && abs(z-rhs.z)<eps;
    }

	double dot(point b) {
	    return (x*b.x + y*b.y + z*b.z);
	}

	void normalize() {
        double div = sqrt(x*x + y*y + z*z);
        if (div != 0) {
            x /= div;
            y /= div;
            z /= div;
        }
    }

    void print() {
        printf("%.7f %.7f %.7f\n", x, y, z);
    }
};

point operator*(double scalar, point& rhs) {
    return rhs * scalar;
}

struct Ray {
    point start;
    point dir;

    Ray() {}

    Ray(point start, point dir) {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};

struct Light {
    point pos;
    double color[3];

    void draw();
};

extern vector<Light*> lights;

struct Object {
    double a, b, c, d, e, f, g, h, i, j;
    point reference_point; // should have x, y, z
    double height, width, length;
    double color[3];
    double coEfficients[4]; // reflection coefficients
    int shine; // exponent term of specular component
    object(){}
    virtual void draw(){}
    void setColor(){}
    void setShine(){}
    void setCoEfficients(){}
    virtual double intersect(Ray* r, double* color, int level) {
        return -1.0;
    }
};

extern vector<Object*> objects;

void printLights() {
    cout << "There are " << lights.size() << " lights." << endl;
    cout << "Also " << objects.size() << " objects." << endl;
}

struct Sphere : public Object {
    int stacks = 50, slices = 50;

    Sphere(point center, double radius) {
        reference_point = center;
        length = radius;
    }

    virtual void draw() {
        point points[stacks+5][slices+5];
        int i,j;
        double h,r;
        //generate points
        for (i=0; i<=stacks; i++) {
            h=length*sin(((double)i/(double)stacks)*(pi/2));
            r=length*cos(((double)i/(double)stacks)*(pi/2));
            for (j=0; j<=slices; j++) {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        glPushMatrix(); {
            glColor3f(color[0], color[1], color[2]);
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            for (i=0; i<stacks; i++) {
                for (j=0; j<slices; j++) {
                    glBegin(GL_QUADS); {
                        //upper half
                        glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                        glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                        glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                        glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                        //lower half
                        glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                        glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                        glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                        glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                    } glEnd();
                }
            }
        } glPopMatrix();
    }

    virtual double intersect(Ray* r, double* color, int level) {
        //coEfficients[0]=1;
        double a = 1.0, b = 2 * r->dir.dot(r->start-reference_point), c = (r->start-reference_point).dot(r->start-reference_point) - length*length;
        double d = b*b - 4*a*c;
        double tMin, tMinus, tPlus;
        if (d < 0) {
            return -1;
        } else {
            d = sqrt(d);
            tMinus = (-b-d)/(2*a);
            tPlus = (-b+d)/(2*a);

            tMin = tMinus <= 0 ? tPlus : tMinus;
        }

        if (level == 0 || tMin <= 0) {
            return tMin;
        } else {
            // update colors
            color[0] = this->color[0]*coEfficients[0];
            color[1] = this->color[1]*coEfficients[0];
            color[2] = this->color[2]*coEfficients[0];

            point p = r->start + tMin*r->dir;
            point n = p - this->reference_point;
            n.normalize();

            for (int i=0; i<lights.size(); i++) {
                Light* l = lights[i];
                Ray r1(l->pos, p-l->pos);
                double t1 = 0, tMin1 = 1000000;
                Object* o;
                Object* oMin = NULL;

                for (int k=0; k<objects.size(); k++) {
                    o = objects[k];
                    t1 = o->intersect(&r1, 0, 0);
                    if (t1 > 0) {
                        if (t1 < tMin1) {
                            tMin1 = t1;
                            oMin = objects[k];
                        }
                    }
                }

                if (oMin == this) {
                    double dist = sqrt((p-r1.start).dot(p-r1.start));
                    //point p1 = r1.start + tMin1*r1.dir;
                    if (abs(dist-tMin1) < eps) {
                        // light hits the object at the intersection point unobstructed
                        r1.dir = -r1.dir;
                        // diffuse
                        double lambert = r1.dir.dot(n);
                        lambert = lambert < 0 ? 0 : lambert;
                        color[0] += this->color[0]*coEfficients[1]*l->color[0]*lambert;
                        color[1] += this->color[1]*coEfficients[1]*l->color[1]*lambert;
                        color[2] += this->color[2]*coEfficients[1]*l->color[2]*lambert;

                        //specular
                        point rr = 2*r1.dir.dot(n)*n - r1.dir;
                        rr.normalize();
                        double phong = rr.dot(-r->dir);
                        phong = phong < 0 ? 0 : phong;
                        color[0] += coEfficients[2]*l->color[0]*pow(phong, shine);
                        color[1] += coEfficients[2]*l->color[1]*pow(phong, shine);
                        color[2] += coEfficients[2]*l->color[2]*pow(phong, shine);
                    }
                }
            }

            if (level < recursion_level) {
                Ray refR;
                refR.start = p;
                refR.dir = r->dir - 2*r->dir.dot(n)*n;
                refR.dir.normalize();
                refR.start = refR.start + 0.001*refR.dir;

                double* refColor = new double[3];
                double tMinR = 1000000, tR = -1;
                Object* oMinR = 0;
                Object* oR;
                for (int k=0; k<objects.size(); k++) {
                    oR = objects[k];
                    tR = oR->intersect(&refR, refColor, 0);
                    if (tR > 0) {
                        if (tR < tMinR) {
                            tMinR = tR;
                            oMinR = objects[k];
                        }
                    }
                }

                if (oMinR) {
                    tR = oMinR->intersect(&refR, refColor, level+1);
                    color[0] += refColor[0] * coEfficients[3];
                    color[1] += refColor[1] * coEfficients[3];
                    color[2] += refColor[2] * coEfficients[3];
                }

                delete[] refColor;
            }

            return tMin;
        }
    }
};

void Light::draw() {
    Sphere s(pos, 2);
    s.color[0] = color[0];
    s.color[1] = color[1];
    s.color[2] = color[2];
    s.draw();
}

struct Triangle : public Object {
    point p1, p2, p3;
    Triangle(point p1, point p2, point p3) {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }

    double det(point a, point b) {
        return a.x*b.y - a.y*b.x;
    }

    double det3(point c1, point c2, point c3) {
        // 3 columns
        return c1.x*(c2.y*c3.z-c3.y*c2.z) + c2.x*(c3.y*c1.z-c1.y*c3.z) + c3.x*(c1.y*c2.z-c2.y*c1.z);
    }

    virtual void draw() {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES); {
            glVertex3f(p1.x, p1.y, p1.z);
            glVertex3f(p2.x, p2.y, p2.z);
            glVertex3f(p3.x, p3.y, p3.z);
        } glEnd();
    }

    virtual double intersect(Ray* r, double* color, int level) {
        point p4 = p2-p1, p5 = p3-p1;
        point n = p4*p5;
        n.normalize();

        point c1(p1-p2), c2(p1-p3), c3(r->dir), c4(p1-r->start);
        double A = det3(c1, c2, c3);
        if (A == 0) return -1;

        double beta = det3(c4, c2, c3) / A;
        double gamma = det3(c1, c4, c3) / A;
        double t = det3(c1, c2, c4) / A;

        if (!((beta+gamma)<1 && beta>0 && gamma>0)) {
            t = -1;
        }

        if (level == 0) {
            return t;
        } else {
            if (t > 0) {
                point p = r->start + t*r->dir;
                color[0] = this->color[0]*coEfficients[0];
                color[1] = this->color[1]*coEfficients[0];
                color[2] = this->color[2]*coEfficients[0];

                for (int i=0; i<lights.size(); i++) {
                    Light* l = lights[i];
                    Ray r1(l->pos, p-l->pos);
                    double t1 = 0, tMin1 = 1000000;
                    Object* o;
                    Object* oMin=NULL;

                    for (int k=0; k<objects.size(); k++) {
                        o = objects[k];
                        t1 = o->intersect(&r1, 0, 0);
                        if (t1 > 0) {
                            if (t1 < tMin1) {
                                tMin1 = t1;
                                oMin = objects[k];
                            }
                        }
                    }

                    if (oMin == this) {
                        point p1 = r1.start + tMin1*r1.dir;

                        if (p==p1) {
                            // light hits the object at the intersection point unobstructed
                            r1.dir = -r1.dir;
                            // diffuse
                            double lambert = r1.dir.dot(n);
                            lambert = lambert < 0 ? 0 : lambert;
                            color[0] += this->color[0]*coEfficients[1]*l->color[0]*lambert;
                            color[1] += this->color[1]*coEfficients[1]*l->color[1]*lambert;
                            color[2] += this->color[2]*coEfficients[1]*l->color[2]*lambert;

                            //specular
                            point rr = 2*r1.dir.dot(n)*n - r1.dir;
                            rr.normalize();
                            double phong = rr.dot(-r->dir);
                            phong = phong < 0 ? 0 : phong;
                            color[0] += coEfficients[2]*l->color[0]*pow(phong, shine);
                            color[1] += coEfficients[2]*l->color[1]*pow(phong, shine);
                            color[2] += coEfficients[2]*l->color[2]*pow(phong, shine);
                        }
                    }
                }

                if (level < recursion_level) {
                    Ray refR;
                    refR.start = p;
                    refR.dir = r->dir - 2*r->dir.dot(n)*n;
                    refR.dir.normalize();
                    refR.start = refR.start + 0.001*refR.dir;

                    double* refColor = new double[3];
                    double tMinR = 1000000, tR = -1;
                    Object* oMinR = 0;
                    Object* oR;
                    for (int k=0; k<objects.size(); k++) {
                        oR = objects[k];
                        tR = oR->intersect(&refR, refColor, 0);
                        if (tR > 0) {
                            if (tR < tMinR) {
                                tMinR = tR;
                                oMinR = objects[k];
                            }
                        }
                    }

                    if (oMinR) {
                        tR = oMinR->intersect(&refR, refColor, level+1);
                        color[0] += refColor[0] * coEfficients[3];
                        color[1] += refColor[1] * coEfficients[3];
                        color[2] += refColor[2] * coEfficients[3];
                    }

                    delete[] refColor;
                }
            }

            return t;
        }
    }
};

struct Floor: public Object {
    bool show_texture = false;
    double floor_width;
    bitmap_image* image;
    Floor(double floor_width, double tile_width) {
        this->floor_width = floor_width;
        reference_point = {-floor_width/2.0, -floor_width/2.0, 0};
        length = tile_width;

        image = new bitmap_image("H:\\Academic\\4-1\\CSE 409 (Graphics)\\Lab\\Offline 3\\Ray Tracing\\Ray Tracing\\texture.bmp");
    }

    virtual void draw() {
        int cnt = ceil(floor_width/length);
        int tile_color = 1;
        for (int i=0; i<cnt; i++) {
            for (int j=0; j<cnt; j++) {
                tile_color = (i+j)%2;
                glColor3f(tile_color, tile_color, tile_color);
                glBegin(GL_QUADS); {
                    point r = {reference_point.x+j*length, reference_point.y+i*length, reference_point.z};
                    glVertex3f(r.x, r.y, r.z);
                    glVertex3f(r.x+length, r.y, r.z);
                    glVertex3f(r.x+length, r.y+length, r.z);
                    glVertex3f(r.x, r.y+length, r.z);
                } glEnd();
            }
        }
    }

    virtual double intersect(Ray* r, double* color, int level) {
        unsigned char rp, gp, bp;
        point n(0, 0, 1);
        double t = -(n.dot(r->start))/(n.dot(r->dir));
        if (t <= 0) {
            return t;
        } else {
            point p;
            p = r->start + r->dir*t;
            int i = (p.y-reference_point.y)/length;
            int j = (p.x-reference_point.x)/length;
            int tileColor = (i+j) % 2;
            int cnt = ceil(floor_width/length);
            if (i>=0 && i<cnt && j>=0 && j<cnt) {
                if (level == 0) return t;

                if (tileColor || !show_texture) {
                    color[0] = color[1] = color[2] = tileColor * this->coEfficients[0];
                } else {
                    unsigned int image_i = (p.y-reference_point.y), image_j = (p.x-reference_point.x);
                    image_i %= 20;
                    image_j %= 20;
                    image->get_pixel(image_i, image_j, rp, gp, bp);
                    color[0] = (rp*1.0/255.0) * this->coEfficients[0];
                    color[1] = (gp*1.0/255.0) * this->coEfficients[0];
                    color[2] = (bp*1.0/255.0) * this->coEfficients[0];
                }

            } else {
                return -1;
            }

            for (int ii=0; ii<lights.size(); ii++) {
                Light* l = lights[ii];
                Ray r1(l->pos, p-l->pos);
                double t1 = 0, tMin1 = 1000000;
                Object* o;
                Object* oMin=NULL;

                for (int k=0; k<objects.size(); k++) {
                    o = objects[k];
                    t1 = o->intersect(&r1, 0, 0);
                    if (t1 > 0) {
                        if (t1 < tMin1) {
                            tMin1 = t1;
                            oMin = objects[k];
                        }
                    }
                }

                if (oMin == this) {
                    point p1 = r1.start + tMin1*r1.dir;

                    if (p==p1) {
                        // light hits the object at the intersection point unobstructed
                        r1.dir = -r1.dir;
                        // diffuse
                        double lambert = r1.dir.dot(n);
                        lambert = lambert < 0 ? 0 : lambert;
                        color[0] += ((tileColor || !show_texture) ? tileColor : (rp*1.0/255))*coEfficients[1]*l->color[0]*lambert;
                        color[1] += ((tileColor || !show_texture) ? tileColor : (gp*1.0/255))*coEfficients[1]*l->color[1]*lambert;
                        color[2] += ((tileColor || !show_texture) ? tileColor : (bp*1.0/255))*coEfficients[1]*l->color[2]*lambert;

                        //specular
                        point rr = 2*r1.dir.dot(n)*n - r1.dir;
                        rr.normalize();
                        double phong = rr.dot(-r->dir);
                        phong = phong < 0 ? 0 : phong;
                        color[0] += coEfficients[2]*l->color[0]*pow(phong, shine);
                        color[1] += coEfficients[2]*l->color[1]*pow(phong, shine);
                        color[2] += coEfficients[2]*l->color[2]*pow(phong, shine);
                    }
                }
            }

            if (level < recursion_level) {
                Ray refR;
                refR.start = p;
                refR.dir = r->dir - 2*r->dir.dot(n)*n;
                refR.dir.normalize();
                refR.start = refR.start + 0.001*refR.dir;

                double* refColor = new double[3];
                double tMinR = 1000000, tR = -1;
                Object* oMinR = 0;
                Object* oR;
                for (int k=0; k<objects.size(); k++) {
                    oR = objects[k];
                    tR = oR->intersect(&refR, refColor, 0);
                    if (tR > 0) {
                        if (tR < tMinR) {
                            tMinR = tR;
                            oMinR = objects[k];
                        }
                    }
                }

                if (oMinR) {
                    tR = oMinR->intersect(&refR, refColor, level+1);
                    color[0] += refColor[0] * coEfficients[3];
                    color[1] += refColor[1] * coEfficients[3];
                    color[2] += refColor[2] * coEfficients[3];
                }

                delete[] refColor;
            }

            return t;
        }
    }
};

struct General : public Object {
    bool eligible(Ray* r, double t) {
        if (t < 0) return false;

        point p = r->start + t*r->dir;
        if (length != 0) {
            if (abs(p.x-reference_point.x) > length) {
                return false;
            }
        }
        if (width != 0) {
            if (abs(p.y-reference_point.y) > width) {
                return false;
            }
        }
        if (height != 0) {
            if (abs(p.z-reference_point.z) > height) {
                return false;
            }
        }

        return true;
    }

    virtual double intersect(Ray* r, double* color, int level) {
        //coEfficients[0]=1;
        double aq = a*r->dir.x*r->dir.x + b*r->dir.y*r->dir.y + c*r->dir.z*r->dir.z + d*r->dir.x*r->dir.y
                    + f*r->dir.y*r->dir.z + e*r->dir.z*r->dir.x;
        double bq = 2*a*r->start.x*r->dir.x + 2*b*r->start.y*r->dir.y + 2*c*r->start.z*r->dir.z +
                    d*(r->start.x*r->dir.y + r->start.y*r->dir.x) + f*(r->start.y*r->dir.z + r->start.z*r->dir.y)
                    + e*(r->start.z*r->dir.x + r->start.x*r->dir.z) + g*r->dir.x + h*r->dir.y + i*r->dir.z;
        double cq = a*r->start.x*r->start.x + b*r->start.y*r->start.y + c*r->start.z*r->start.z + d*r->start.x*r->start.y
                    + f*r->start.y*r->start.z + e*r->start.z*r->start.x + g*r->start.x + h*r->start.y + i*r->start.z + j;

        double tMinus, tPlus, tMin;
        if (aq == 0) {
            tMin = -cq/bq;
            if (tMin <= 0) {
                return tMin;
            }
        } else {
            double dq = bq*bq - 4*aq*cq;
            if (dq < 0) {
                return -1;
            } else {
                dq = sqrt(dq);
                tMinus = (-bq-dq)/(2*aq);
                tPlus = (-bq+dq)/(2*aq);

                if (!eligible(r, tMinus) && !eligible(r, tPlus)) {
                    return -1;
                }
                tMin = !eligible(r, tMinus) ? tPlus : tMinus;
            }
        }

        if (level == 0) {
            return tMin;
        } else {
            if (tMin >= 0) {
                color[0] = this->color[0]*coEfficients[0];
                color[1] = this->color[1]*coEfficients[0];
                color[2] = this->color[2]*coEfficients[0];

                point p = r->start + tMin*r->dir;
                point n;
                n.x = 2*a*p.x + d*p.y + e*p.z + g;
                n.y = 2*b*p.y + d*p.x + f*p.z + h;
                n.z = 2*c*p.z + e*p.x + f*p.y + i;
                n.normalize();

                for (int ii=0; ii<lights.size(); ii++) {
                    Light* l = lights[ii];
                    Ray r1(l->pos, p-l->pos);
                    double t1 = 0, tMin1 = 1000000;
                    Object* o;
                    Object* oMin=NULL;

                    for (int k=0; k<objects.size(); k++) {
                        o = objects[k];
                        t1 = o->intersect(&r1, 0, 0);
                        if (t1 > 0) {
                            if (t1 < tMin1) {
                                tMin1 = t1;
                                oMin = objects[k];
                            }
                        }
                    }

                    if (oMin == this) {
                        point p1 = r1.start + tMin1*r1.dir;

                        if (p==p1) {
                            // light hits the object at the intersection point unobstructed
                            r1.dir = -r1.dir;
                            // diffuse
                            double lambert = r1.dir.dot(n);
                            lambert = lambert < 0 ? 0 : lambert;
                            color[0] += this->color[0]*coEfficients[1]*l->color[0]*lambert;
                            color[1] += this->color[1]*coEfficients[1]*l->color[1]*lambert;
                            color[2] += this->color[2]*coEfficients[1]*l->color[2]*lambert;

                            //specular
                            point rr = 2*r1.dir.dot(n)*n - r1.dir;
                            rr.normalize();
                            double phong = rr.dot(-r->dir);
                            phong = phong < 0 ? 0 : phong;
                            color[0] += coEfficients[2]*l->color[0]*pow(phong, shine);
                            color[1] += coEfficients[2]*l->color[1]*pow(phong, shine);
                            color[2] += coEfficients[2]*l->color[2]*pow(phong, shine);
                        }
                    }
                }

                if (level < recursion_level) {
                    Ray refR;
                    refR.start = p;
                    refR.dir = r->dir - 2*r->dir.dot(n)*n;
                    refR.dir.normalize();
                    refR.start = refR.start + 0.001*refR.dir;

                    double* refColor = new double[3];
                    double tMinR = 1000000, tR = -1;
                    Object* oMinR = 0;
                    Object* oR;
                    for (int k=0; k<objects.size(); k++) {
                        oR = objects[k];
                        tR = oR->intersect(&refR, refColor, 0);
                        if (tR > 0) {
                            if (tR < tMinR) {
                                tMinR = tR;
                                oMinR = objects[k];
                            }
                        }
                    }

                    if (oMinR) {
                        tR = oMinR->intersect(&refR, refColor, level+1);
                        color[0] += refColor[0] * coEfficients[3];
                        color[1] += refColor[1] * coEfficients[3];
                        color[2] += refColor[2] * coEfficients[3];
                    }

                    delete[] refColor;
                }
            }

            return tMin;
        }
    }
};


