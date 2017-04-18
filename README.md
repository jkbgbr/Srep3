# SRep

SRep - ShapeRepresentation - is aimed to be an element of a larger program.
The problem to be solved is identifying closed volumes defined by cylindrical shells and heads. The volumes may be
neighbours, may be fully or partially encapsulated in each-other; see the tests for examples.
The way of solving the problem is simplifying the 3D problem to 2D by taking a section, and finding the relation
between the polygons representing the original volumes.
The aim of the package is to provide abstact tools that can reliably identify neighbouring elements,
internal points etc. with using as little external dependencies as possible.
Tests are provided as examples and to make sure the results are correct.

The method is:
- polygons are created either generated from bottom-up (Points make up LineSegments make up PolygonsLines make up
PolygonGroups) or using some elements from the module SR_vessel. These represent the 3D shells and together they make up (span) closed areas - polygons.
- sets of non-overlapping polygons is found by triangulating the defined closed areas
- returned is a dict describing the relation of the closed areas and the polygons spanning these.

There are some methods used that may be interesting for other use cases:
- straight line segment relative positions
- convex hull of 2D point cloud
- triangulation algorithms (delaunay, ear clipping)

Note: the algorithms and methods used here are by far not the most efficient ones, but they get the job done.

**Coordinate system**
- origin in the lower left corner
- x axis horizontal, growing eastwards
- y axis vertical, growing northwards
- angle is measured in degrees from the x axis, CCW is positive

**Point**

Point objects are considered 0D Entities, defined by a tuple. Coordinates defined as int or float are accepted, converted to float.

`p1 = Point((1, 2)) == Point((1.0, 2.0))`

The SR_point module provides basic function to validate/modify the Point objects.

**LineSegment**

LineSegment objects are directed 1D Entities, defined by two endpoints, each a valid Point object.
LineSegments are considered as a section between the two endpoints on an otherwise infinite line.
The direction of the Segments is from endpoint "i" to "j" with i, j being the first and second (valid) Points
used to define the LineSegment.
The interior of a LineSegment is the line itself including the endpoints.

`p1 = Point((1.0, 2.0))`

`p2 = Point((1, 5))`

`ls = LineSegment(p1, p2)`

A set of methods is provided to evaluate the relative positions of segment lines (angle, parallelity, collinearity,
intersection point, overlap etc.) and lines to points (distance, sidedness of point relative to line etc).
This is by far not a full set of 2D tools as in a CAD program.

**IntexPoint**

IntexPoint objects are Points that have a parent LineSegment to which they belong, but are themselves no part of that LineSegment.

`ip1 = IntexPoint((1, 5), ls)`

**PolygonLine**

PolygonLine objects are sets of LineSegments that form an open or closed polygon.

`p1 = Point((1.0, 2.0))`

`p2 = Point((1, 5))`

`p3 = Point((4, 6))`

`ls1 = LineSegment(p1, p2)`

`ls2 = LineSegment(p2, p3)`

`pl1 = PolygonLine(ls1, ls2)`

**Triangle**

Triangle objects are three-sided PolygonLine objects.

A set of methods is provided to evaluate the relative position of triangles (congruence, overlapping, touching etc.), and a smaller
set of methods for arbitrary polygons.

The SR_polygon module provides functions to create/modify basic shapes (square, arc) and modify them (translation, rotation).
It also enables to define IntexPoints as known internal/external Points to the segments making up the PolygonLine.
These are used in evaluating how the neighbouring areas relate to each other.

**PolygonGroup**

These are made up by PolygonLines. These are used to represent the 3D objects in 2D and find the volumes.

**Graph, BehalterGraph**

As the polygons may consist of lots of segments to provide a realistic approximation of the real geometry and some algorithms
may be slow, a Graph is a simplyfied (or equal) representation of the 2D geometry. Finding the areas and their relations is done
using the graph leading to acceptable running times.

**SR_hauptteile module**

This module contains the high-level classes to be used to define the polygongroups to solve the problems this whole package was created for.