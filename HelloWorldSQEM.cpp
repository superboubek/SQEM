/*
 *     SQEM - Spherical Quadric Error Metric
 *
 *     Authors:  Jean-Marc THIERY, Emilie GUY and Tamy BOUBEKEUR
 *
 *     Copyright © 2012-2013 Telecom ParisTech - CNRS LTCI - Institut Mines Telecom
 *              All rights reserved
 *
 * This file is a part of the standalone implementation of the 
 * 
 * Spherical Quadric Error Metric.
 * 
 * For more information or if you use this file and need to reference it, please
 * 
 * refer to the following publication:
 *
 *    Sphere-Meshes: Shape Approximation using Spherical Quadric Error Metrics
 *    
 *    Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur
 *    
 *    ACM Transaction on Graphics (Proc. SIGGRAPH Asia 2013), 32(6), Art. 178
 *
 *    http://www.telecom-paristech.fr/~boubek/papers/SphereMeshes/
 * 
 * SQEM is free software: you can redistribute it and/or modify
 *
 * it under the terms of the GNU Lesser General Public License as published by
 *
 * the Free Software Foundation, either version 3 of the License, or
 *
 * (at your option) any later version.
 *
 * SQEM is distributed in the hope that it will be useful,
 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cstdlib>
#include "SQEM.h"

using namespace std;

class Point {
public:
	inline Point (float x = 0.f, float y = 0.f, float z = 0.f) {
		p[0] = x;
		p[1] = y;
		p[2] = z;
	}
	inline float & operator[] (unsigned int i) { return p[i]; }
	inline const float & operator[] (unsigned int i) const { return p[i]; }
	inline Point& operator= (const Point & P) {
	  p[0] = P[0];
	  p[1] = P[1];
	  p[2] = P[2];
	  return (*this);
	};

	inline Point& operator+= (const Point & P) {
	  p[0] += P[0];
	  p[1] += P[1];
	  p[2] += P[2];
	  return (*this);
	};

	inline Point& operator-= (const Point & P) {
	  p[0] -= P[0];
	  p[1] -= P[1];
	  p[2] -= P[2];
	  return (*this);
	};

	inline Point& operator*= (const Point & P) {
	  p[0] *= P[0];
	  p[1] *= P[1];
	  p[2] *= P[2];
	  return (*this);
	};

	inline Point& operator*= (float s) {
	  p[0] *= s;
	  p[1] *= s;
	  p[2] *= s;
	  return (*this);
	};

	inline Point& operator/= (const Point & P) {
	  p[0] /= P[0];
	  p[1] /= P[1];
	  p[2] /= P[2];
	  return (*this);
	};

	inline Point& operator/= (float s) {
	  p[0] /= s;
	  p[1] /= s;
	  p[2] /= s;
	  return (*this);
	};

	inline Point operator+ (const Point & P) const {
	  Point res;
	  res[0] = p[0] + P[0];
	  res[1] = p[1] + P[1];
	  res[2] = p[2] + P[2];
	  return (res); 
	};

	inline Point operator- (const Point & P) const {
	  Point res;
	  res[0] = p[0] - P[0];
	  res[1] = p[1] - P[1];
	  res[2] = p[2] - P[2];
	  return (res); 
	};

	inline Point operator- () const {
	  Point res;
	  res[0] = -p[0];
	  res[1] = -p[1];
	  res[2] = -p[2];
	  return (res); 
	};

	inline Point operator* (const Point & P) const {
	  Point res;
	  res[0] = p[0] * P[0];
	  res[1] = p[1] * P[1];
	  res[2] = p[2] * P[2];
	  return (res); 
	};

	inline Point operator* (float s) const {
	  Point res;
	  res[0] = p[0] * s;
	  res[1] = p[1] * s;
	  res[2] = p[2] * s;
	  return (res); 
	};

	inline Point operator/ (const Point & P) const {
	  Point res;
	  res[0] = p[0] / P[0];
	  res[1] = p[1] / P[1];
	  res[2] = p[2] / P[2];
	  return (res); 
	};

	inline Point operator/ (float s) const {
	  Point res;
	  res[0] = p[0] / s;
	  res[1] = p[1] / s;
	  res[2] = p[2] / s;
	  return (res); 
	};
	
	inline bool operator == (const Point & a) const {
	  return(p[0] == a[0] && p[1] == a[1] && p[2] == a[2]);
	};
	
	inline bool operator != (const Point & a) const {
	  return(p[0] != a[0] || p[1] != a[1] || p[2] != a[2]);
	};
		
	inline bool operator < (const Point & a) const {
	  return(p[0] < a[0] && p[1] < a[1] && p[2] < a[2]);
	};
		
	inline bool operator >= (const Point &a) const {
	  return(p[0] >= a[0] && p[1] >= a[1] && p[2] >= a[2]);
	};
private:
	float p[3];
};

inline Point operator * (float s, const Point & p) {
	return (p * s);
}

inline float dot (const Point & a, const Point & b) {
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline float length (const Point & x) {
	return sqrt (dot (x, x));
}

inline Point normalize (const Point & x) {
    Point n (x);
    n /= length (n);
    return n;
}

int main (int argc, char ** argv) {
	cout << "-----------------------------" << endl;
	cout << "SQEM HelloWorld Test Program." << endl;
	cout << "  Step 1: initializing a SQEM with the six planes of a unit cube." << endl;
	SQEM sqem (Point (1.f, 0.f, 0.f), Point (1.f, 0.f, 0.f));
	sqem += SQEM (Point (-1.f, 0.f, 0.f), Point (-1.f, 0.f, 0.f));
	sqem += SQEM (Point (0.f, 1.f, 0.f), Point (0.f, 1.f, 0.f));
	sqem += SQEM (Point (0.f, -1.f, 0.f), Point (0.f, -1.f, 0.f));
	sqem += SQEM (Point (0.f, 0.f, 1.f), Point (0.f, 0.f, 1.f));
	sqem += SQEM (Point (0.f, 0.f, -1.f), Point (0.f, 0.f, -1.f));
	Point sphereCenter;
	float sphereRadius;
	const float MAX_RADIUS = 10.f;
	cout << "  Step 2: optimizing a sphere which minimze the SQEM." << endl;
	sqem.minimize (sphereCenter, sphereRadius, Point (-MAX_RADIUS, -MAX_RADIUS, -MAX_RADIUS), Point (MAX_RADIUS, MAX_RADIUS, MAX_RADIUS));
	cout << "    Result: optimal sphere centered at [" << sphereCenter[0] << ", " << sphereCenter[1] << ", " << sphereCenter[3] << "], with radius " << sphereRadius << "." << endl;
	cout << "  Step 3: reinitializing the SQEM using the six planes of a smaller shifted cube." << endl;
	sqem.setFromPlan (Point (2.5f, 2.f, 0.f), Point (1.f, 0.f, 0.f));
    sqem += SQEM (Point (1.5f, 2.f, 0.f), Point (-1.f, 0.f, 0.f));
	sqem += SQEM (Point (2.f, 2.5f, 0.f), Point (0.f, 1.f, 0.f));
	sqem += SQEM (Point (2.f, 1.5f, 0.f), Point (0.f, -1.f, 0.f));
	sqem += SQEM (Point (2.f, 2.f, 0.5f), Point (0.f, 0.f, 1.f));
	sqem += SQEM (Point (2.f, 2.f, -0.5f), Point (0.f, 0.f, -1.f));
	cout << "  Step 4: optimizing a sphere which minimze the SQEM." << endl;
	sqem.minimize (sphereCenter, sphereRadius, Point (2.f-MAX_RADIUS, 2.f-MAX_RADIUS, 0.f), Point (2.f+MAX_RADIUS, 2.f+MAX_RADIUS, 0.f));
	cout << "    Result: optimal sphere centered at [" << sphereCenter[0] << ", " << sphereCenter[1] << ", " << sphereCenter[3] << "], with radius " << sphereRadius << "." << endl;
	cout << "  Step 5: reinitializing the SQEM using 1000 planes randomly sampled on the unit sphere." << endl;
	sqem.setFromPlan (Point (1.0f, 0.f, 0.f), Point (1.f, 0.f, 0.f));
	srand (0);
	for (unsigned int i = 0; i < 1000; i++) {
		Point r (float (std::rand () - (RAND_MAX/2)), float (std::rand () - (RAND_MAX/2)), float (std::rand () - (RAND_MAX/2)));
		Point p = normalize (r);
		sqem += SQEM (p, p);
	}
	cout << "  Step 6: optimizing a sphere which minimze the SQEM." << endl;
	sqem.minimize (sphereCenter, sphereRadius, Point (-MAX_RADIUS, -MAX_RADIUS, -MAX_RADIUS), Point (MAX_RADIUS, MAX_RADIUS, MAX_RADIUS));
	cout << "    Result: optimal sphere centered at [" << sphereCenter[0] << ", " << sphereCenter[1] << ", " << sphereCenter[3] << "], with radius " << sphereRadius << "." << endl;
	cout << "-----------------------------" << endl;
	return 0;
}