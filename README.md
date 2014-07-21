## Seismic image processing for geologic faults

This repository contains computer programs written and used by 
[Dave Hale](http://inside.mines.edu/~dhale) 
to create examples shown in his SEG/AAPG Distinguished Lecture in Fall 2014,
entitled 
[3D seismic image processing for interpretation of faults and horizons]
(http://www.seg.org/education/lectures-courses/distinguished-lecturers/fall2014/hale-abstract).

This software depends on that in the [Mines Java Toolkit
(JTK)](https://github.com/dhale/jtk/). If you want to do more than browse the
source code, you must first download and build the Mines JTK using
[Gradle](http://www.gradle.org). The build process for software in
this repository is the same.

Like the Mines JTK, this is a toolkit for computer programmers. It is not a
complete system for seismic interpretation of geologic faults. Others
(including commercial software companies) have built such systems using
earlier versions of one or more of the software tools provided in this
repository.

Although based on [research first published in 2012]
(http://inside.mines.edu/~dhale/research.html)
most of the software in this repository was written in the summer of 2014
as Dave prepared for the upcoming lecture tour. This software will change
during the course of that tour in response to questions and discussions
inspired by the lecture.

### Summary

Here are brief descriptions of key components:

#### FaultScanner
Creates images of fault likelihood, strike and dip. This software also
includes methods for thinning such images and for using them to smooth seismic
images along reflectors, while enhancing discontinuities at faults. An earlier
version of the computationally efficient scan over fault strikes and dips has
been used in both commercial and in-house systems for seismic interpretation.

#### FaultSkinner (new in 2014)
Constructs fault skins comprised of fault cells. A fault cell is a point
located on a fault, with an associated fault likelihood, strike and dip. A
fault skin is a linked collection of fault cells used to walk within a seismic
image, up and down along fault curves (tangent to fault dip) and left and
right along fault traces (tangent to fault strike). In contrast to surface
meshes often used to represent faults, a fault skin is a lightweight data
structure designed specifically to facilitate image processing.

#### FaultSlipper (new in 2014)
Estimates fault dip slips using seismic amplitudes on opposite sides of fault
skins. (Does not estimate strike slips.) Also includes methods to construct a
complete field of fault dip-slip vectors, and to use those vectors to undo
faulting apparent in a seismic image. *Unfaulting* a seismic image is one way
to assess the accuracy of the estimated dip slips.

#### FakeData (new in 2014)
Generates synthetic seismic images useful for research. In testing new ideas,
it often helps to know the answers.

---
Copyright (c) 2014, Dave Hale. All rights reserved.
This software and accompanying materials are made available under the terms of
the [Common Public License - v1.0](http://www.eclipse.org/legal/cpl-v10.html),
which accompanies this distribution.
