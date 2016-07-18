# SQEM
Standalone Spherical Quadric Error Metric C++ implementation, introduced with Sphere-Meshes (Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur, ACM SIGGRAPH Asia 2013).

## VERSION
Version 1.0

## DESCRIPTION
The SQEM class allows to compute a sphere which approximates a set of planes in 3D.

A Spherical Quadric Error Metric (or SQEM) models the squared distance between a sphere and plane 
or a set of plane. Its minimizer is a sphere that approximates the set of plane. 
See Sphere-Meshes: Shape Approximation using Spherical Quadric Error Metrics, 
by Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur 
(ACM Transaction on Graphics (Proc. SIGGRAPH Asia 2013), 32(6), Art. 178) 
for more information.

The SQEM implementation is provided as a single header for easy insertion in your projects. 
It has been tested under Linux Ubuntu 14.04 and Windows 7, in both case using GCC.

## INSTALL

To compile the demonstration program in a terminal:
    g++ -O2 -o HelloWorldSQEM HelloWorldSQEM.cpp
    ./HelloWorldSQEM

## CODE DOCUMENTATION

Just run Doxygen (www.doxygen.org) on the Doxyfile of SQEM to generate an html
documentation of the source code, which will be located in the doc
directory:
 
  doxygen Doxyfile

## CITATION

Please cite the following paper when using SimSelect:

Sphere-Meshes: Shape Approximation using Spherical Quadric Error Metrics
Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur
ACM Transaction on Graphics (Proc. SIGGRAPH Asia 2013), 32(6), Art. 178

## CREDITS and LICENCE

     SQEM - Spherical Quadric Error Metric

     Authors:  Jean-Marc THIERY, Emilie GUY and Tamy BOUBEKEUR

     Copyright © 2012-2013 LTCI, CNRS, Telecom ParisTech, University Paris-Saclay
              All rights reserved

 This file is a part of the standalone implementation of the 
 
 Spherical Quadric Error Metric.
 
 For more information or if you use this file and need to reference it, please
 
 refer to the following publication:

    Sphere-Meshes: Shape Approximation using Spherical Quadric Error Metrics
    
    Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur
    
    ACM Transaction on Graphics (Proc. SIGGRAPH Asia 2013), 32(6), Art. 178

    http://www.telecom-paristech.fr/~boubek/papers/SphereMeshes/
 
 SQEM is free software: you can redistribute it and/or modify

 it under the terms of the GNU Lesser General Public License as published by

 the Free Software Foundation, either version 3 of the License, or

 (at your option) any later version.

 SQEM is distributed in the hope that it will be useful,

 but WITHOUT ANY WARRANTY; without even the implied warranty of

 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the

 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License

 along with this program. If not, see <http://www.gnu.org/licenses/>.

## CONTACT

Webpage: http://www.telecom-paristech.fr/~boubek/papers/SphereMeshes
Contact: jean-marc.thiery@telecom-paristech.fr, emilie.guy@telecom-paristech.fr, tamy.boubekeur@telecom-paristech.fr

