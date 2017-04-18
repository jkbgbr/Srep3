# -*- coding: utf-8 -*-

from geometry.SR_line import LineSegment
from geometry.SR_point import Point, IntexPoint
from geometry.SR_polygon import PolygonLine, square, convex_hull, Triangle, remove_collinear_points
import unittest
import math
import pprint as pp
from geometry import PURGE
import itertools

EPS = 1e-10  # accepted numerical error


class TestPolygon(unittest.TestCase):

    def setUp(self):
        self.l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.l2 = LineSegment(Point((0, 1)), Point((2, 3)))
        self.l3 = LineSegment(Point((4, 4)), Point((5, 5)))

    def test_definition(self):
        self.assertRaises(Exception, PolygonLine.__init__)

    def test_loose_subset(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        l2 = LineSegment(Point((0, 1)), Point((2, 3)))
        l3 = LineSegment(Point((2, 3)), Point((5, 5)))
        ss1 = LineSegment(Point((0, 0)), Point((0, 1)))
        ss2 = LineSegment(Point((2, 3)), Point((0, 1)))
        ss3 = LineSegment(Point((20, 30)), Point((0, 10)))
        pl1 = PolygonLine([l1, l2, l3])
        pl2 = PolygonLine([ss1, ss2])
        pl3 = PolygonLine([ss1, ss3])
        self.assertEqual(pl1.loose_subset(pl2), True)
        self.assertEqual(pl1.loose_subset(pl3), False)

    def test_magic_functions(self):
        # testing the magic functions
        l1 = LineSegment()
        self.assertNotEqual(l1, LineSegment())  # these are both non-valid
        l3 = LineSegment(Point((0, 1)), Point((2, 4)))
        pl1 = PolygonLine([self.l1, self.l2])  # pl1
        pl2 = PolygonLine([self.l1, self.l2])  # same as pl1
        pl3 = PolygonLine([self.l1, l3])  # not the same as pl1
        self.assertEqual(pl1, pl2)
        self.assertNotEqual(pl1, pl3)
        self.assertEqual(pl1, eval(pl1.__repr__()))

    def test_loose_equality(self):
        # testing loose equality of polygon lines - open polygons only
        p1 = Point((0, 0))
        p2 = Point((1, 1))
        p3 = Point((2, 1.5))
        p4 = Point((3, 1))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p2, p4)
        ls11 = LineSegment(p2, p1)  # reverse von ls1
        pl1 = PolygonLine([ls1, ls2])
        pl2 = PolygonLine([ls1, ls3])
        pl3 = PolygonLine([ls11, ls2])  # first segment is reverse of ls1
        pl4 = PolygonLine([ls2, ls11])  # segment order different, one of them is reversed
        self.assertEqual(pl1 == pl2, False)  # not equal, as a segment is entirely different
        self.assertEqual(pl1 == pl3, False)  # not equal, as one of the segments is reversed
        self.assertEqual(pl1.loose_equal(pl3), True)  # first segment is reversed, loose equal
        self.assertEqual(pl1.loose_equal(pl4), True)  # segment order, segment direction different, loose equal
        pl5 = PolygonLine([ls1])
        self.assertEqual(pl1 == pl5, False)  # not equal lengths

    def test_appending_point_at_end(self):
        # testing equality, adding a Point
        pl1 = PolygonLine([self.l1, self.l2])
        p = Point((2, 5))
        pl1.append_point_at_end(p)
        l3 = LineSegment(Point((2, 3)), p)
        pl2 = PolygonLine([self.l1, self.l2, l3])
        self.assertEqual(pl1, pl2)
        self.assertEqual(pl2, eval(pl2.__repr__()))

    def test_segment_reserved(self):
        # one of the segments
        pl1 = PolygonLine([self.l1, self.l2])
        import copy
        pl2 = copy.deepcopy(pl1)
        pl2._segments[0].reverse()  # last segment reversed
        self.assertNotEqual(pl1, pl2)  # not equal anymore
        self.assertNotEqual(pl2.isvalid, True)  # pl2 not valid

    def test_point_set(self):
        pl1 = PolygonLine([self.l1, self.l2])
        pl1.close()  # adds a third segment
        self.assertTupleEqual(tuple(pl1.point_set), (Point((0, 0)), Point((0, 1)), Point((2, 3))))

    def test_closing(self):
        # creates a 2-segment polygon, closes it, reverses last semgent and checks for continuity,
        # then orders and re-tests
        pl1 = PolygonLine([self.l1, self.l2])
        pl1.close()  # adds a third segment
        self.assertEqual(pl1.isclosed, True)
        self.assertEqual(pl1.iscontinuous, True)
        pl1._segments[-1].reverse()
        self.assertEqual(pl1.iscontinuous, False)
        pl1.order()
        self.assertEqual(pl1.iscontinuous, True)

    def test_validity(self):
        pl1 = PolygonLine([self.l1, self.l2])  # this is OK
        self.assertEqual(pl1.isvalid, True)
        pl1 = PolygonLine([self.l1, self.l3])  # not continuous - missing part
        self.assertEqual(pl1.isvalid, False)
        import copy
        s2 = copy.deepcopy(self.l2)
        s2.reverse()
        pl1 = PolygonLine([self.l1, s2])  # not continuous - 2nd segment reversed
        self.assertEqual(pl1.isvalid, False)
        pl1 = PolygonLine([self.l1])  # just one segment
        self.assertEqual(pl1.isvalid, True)

    def test_open_close(self):
        pl1 = PolygonLine([self.l1, self.l2])
        self.assertEqual(pl1.isclosed, False)
        pl1.close()  # adds a third segment
        self.assertEqual(pl1.isclosed, True)

    def test_closed(self):
        # closed version
        pl1 = PolygonLine([self.l1, self.l2])
        pl2 = pl1.closed_copy()
        self.assertEqual(pl1.isclosed, False)
        self.assertEqual(pl2.isclosed, True)
        self.assertTupleEqual(tuple(pl1.point_set), tuple(pl2.point_set))
        pl3 = pl2.closed_copy()
        self.assertEqual(pl2, pl3)

    def test_rotation_about_origin(self):
        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((2, 3))
        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        pl1 = PolygonLine([l1, l2])
        pl1.close()  # closes
        pl2 = pl1.closed_copy()  # creates a copy of the closed version
        pl2.rotate_about_origin(phi=360)
        from geometry.SR_point import distance
        self.assertTrue(all([distance(x, y) < EPS for x, y in zip(pl1.point_set, pl2.point_set)]))

    def test_rotation(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        l2 = LineSegment(Point((0, 1)), Point((2, 3)))
        pl1 = PolygonLine([l1, l2])
        pl2 = pl1.closed_copy()  # creates a copy and closes it
        pl2.rotate_about_origin(phi=180)  # all coordinates are negated
        self.assertTrue(all([z.x == -w.x for z, w in zip(pl1.point_set, pl2.point_set[::-1])]))
        self.assertTrue(all([z.y == -w.y for z, w in zip(pl1.point_set, pl2.point_set[::-1])]))

    def test_moving(self):
        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((2, 3))
        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        pl1 = PolygonLine([l1, l2])
        pl1.close()  # closes
        pl2 = pl1.closed_copy()  # creates a copy of the closed version
        pl2.move(x0=3, y0=4)  # 3, 4, 5 pythagorean triplet
        from geometry.SR_point import distance
        self.assertTrue(all([distance(x, y) - 5 < EPS for x, y in zip(pl1.point_set, pl2.point_set)]))

    # def test11(self):
    #     # a test for the rotation about arbitrary point
    #     l1 = LineSegment(Point((0, 0)), Point((0, 1)))
    #     l2 = LineSegment(Point((0, 1)), Point((1, 1)))
    #     pl1 = PolygonLine([l1, l2])
    #     pl2 = pl1.closed_copy()  # creates a copy and closes it
    #     pl2.plot()
    #     _center = Point((2, 2))
    #     pl2.rotate_about_xy(x0=_center.x, y0=_center.y, phi=180)  # all coordinates are negated
    #     pl2.plot()
    #     from geometry.SR_point import distance
    #     print([w.x == z.x + 2 * distance(z, _center) for z, w in zip(pl1.point_set, pl2.point_set[::-1])])
    #     self.assertTrue(all([z.x == z.x + 2 * distance(z, _center) for z, w in zip(pl1.point_set, pl2.point_set[::-1])]))
    #     # self.assertTrue(all([z.y == -w.y for z, w in zip(pl1.point_set, pl2.point_set[::-1])]))

    def test_rotation_matrix(self):
        from geometry.SR_polygon import rotationsmatrix
        _phi = 45
        self.assertAlmostEqual(rotationsmatrix(_phi)[0][0], 0.5 * math.sqrt(2), delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[0][1], 0.5 * math.sqrt(2), delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[1][0], - 1 * 0.5 * math.sqrt(2), delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[1][1], 0.5 * math.sqrt(2), delta=EPS)
        _phi = 30
        self.assertAlmostEqual(rotationsmatrix(_phi)[0][0], 0.5 * math.sqrt(3), delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[0][1], 0.5, delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[1][0], - 0.5, delta=EPS)
        self.assertAlmostEqual(rotationsmatrix(_phi)[1][1], 0.5 * math.sqrt(3), delta=EPS)

    def test_purge(self):
        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((2, 3))
        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        pl1 = PolygonLine([l1, l2])
        pl1.purge()
        pl1.purge()
        pl1.close()
        self.assertEqual(pl1.isconvex, True)

        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((0, 1))
        p4 = Point((2, 3))
        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p3, p4)
        pl1 = PolygonLine([l1, l2])
        pl1.purge()
        pl1.purge()
        pl1.close()
        self.assertEqual(pl1.isconvex, True)

    def test_convexity(self):
        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((0, 3))
        p4 = Point((1, 3))
        p5 = Point((1, 0))

        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        l3 = LineSegment(p3, p4)
        l4 = LineSegment(p4, p5)
        pl1 = PolygonLine([l1, l2, l3, l4])
        pl1.close()
        pl1.purge()

        self.assertEqual(pl1.isconvex, True)

    def test_concavity(self):
        p1 = Point((0, 0))
        p2 = Point((0.5, 1))
        p3 = Point((0, 3))
        p4 = Point((1, 3))
        p5 = Point((1, 0))

        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        l3 = LineSegment(p3, p4)
        l4 = LineSegment(p4, p5)
        pl1 = PolygonLine([l1, l2, l3, l4])
        pl1.close()
        pl1.purge()
        self.assertEqual(pl1.isconvex, False)

    def test_finding_segment(self):
        # find segment by point
        p1 = Point((0, 0))
        p2 = Point((0.5, 1))
        p3 = Point((0, 3))
        p4 = Point((1, 3))
        p5 = Point((1, 0))

        l1 = LineSegment(p1, p2)  # defined so that there is no need to purge
        l2 = LineSegment(p2, p3)
        l3 = LineSegment(p3, p4)
        l4 = LineSegment(p4, p5)
        pl1 = PolygonLine([l1, l2, l3, l4])
        self.assertEqual(pl1.find_segment_by_point(point=p2)[0], l2)
        self.assertEqual(pl1.find_segment_by_point(point=Point((12, 15))), [])
        pl1.purge()
        self.assertEqual(pl1.find_segment_by_point(point=p2, whichend='end')[0], l1)
        # todo: case when there is a branch and the point asked for is the branching point

    def test_finding_crossing_lines(self):
        p1 = Point((0, 0))
        p2 = Point((0, 1))
        p3 = Point((1, 1))
        p4 = Point((1, 0))
        l1 = LineSegment(p1, p3)
        l2 = LineSegment(p3, p2)
        l3 = LineSegment(p2, p4)
        l4 = LineSegment(p4, p1)
        pl1 = PolygonLine([l1, l2, l3, l4])
        self.assertEqual(pl1.crossing_lines, [(l1, l3), (l3, l1)])
        self.assertEqual(pl1.isvalid, False)  # should be False b/c crossing lines

    def test_adding_known_point(self):
        # adding internal, external point to convex
        pl1 = square(a=1, x0=0.5, y0=0.5)
        self.assertEqual(pl1.isvalid, True)
        # self.assertEqual(pl1.iscontinuous, True)
        # self.assertEqual(pl1.isconvex, True)
        # pi = Point((0.6, 0.6))
        # pe = Point((0.6, 0.4))
        # # pl1.plot(also_plot=[pi, pe])
        # pl1.add_internal_point(p=pi)
        # self.assertEqual(pi in pl1.internal, True)
        # pl1.add_external_point(p=pe)
        # self.assertEqual(pe in pl1.external, True)
        # self.assertEqual(all([pl1.is_internal_point(p=x) for x in pl1.internal]), True)
        # self.assertEqual(all([pl1.is_external_point(p=x) for x in pl1.external]), True)

    def test_convex_hull_randomized(self):
        import random
        for i in range(1):
            _pointnr = random.randint(10, 100)
            pointset = [(random.uniform(-1, 1) * random.uniform(0, 10) + random.uniform(0, 5), random.uniform(0, 1) * random.uniform(0, 10) + random.uniform(0, 10)) for i in range(_pointnr)]
            ch = convex_hull(pointset)  # creates the random hull

            count = 0
            for p in pointset:
                if not ch.is_internal_point(Point(p)):  # if not internal
                    self.assertEqual(Point(p) in ch.point_set, True)  # then in the convex hull
                else:
                    count += 1

            self.assertEqual(len(ch.point_set) + count, _pointnr)  # no outliers
            self.assertEqual(ch.isvalid, True)
            self.assertEqual(ch.iscontinuous, True)
            self.assertEqual(ch.isconvex, True)
            self.assertEqual(ch.iscw, False)
            self.assertEqual(ch.isccw, True)

        # ch.plot(other_entities=[Point(n) for n in pointset])

    def test_known_point_concave_polygon_ccw(self):
        # adding internal, external point to concave
        pl1 = square(a=1, x0=0.5, y0=0.5)
        pl1._segments[0].i._set_coords(coords=(0.5, 0.7))
        self.assertEqual(pl1.isvalid, True)
        self.assertEqual(pl1.iscontinuous, True)
        self.assertEqual(pl1.isconvex, False)
        self.assertEqual(pl1.iscw, False)
        self.assertEqual(pl1.isccw, True)
        pi = IntexPoint((0.7, 0.7))
        pe = IntexPoint((0.6, 0.4))
        pl1.add_internal_point(p=pi)
        self.assertEqual(pi in pl1.internal, True)
        self.assertEqual(pi in pl1.external, False)
        pl1.add_external_point(p=pe)
        self.assertEqual(pe in pl1.external, True)
        self.assertEqual(pe in pl1.internal, False)
        self.assertEqual(all([pl1.is_internal_point(p=x) for x in pl1.internal]), True)
        self.assertEqual(all([pl1.is_external_point(p=x) for x in pl1.external]), True)

    def test_known_point_concave_polygon_cw(self):
        pl1 = square(a=1, x0=0.5, y0=0.5, cw=True)
        pl1._segments[0].i._set_coords(coords=(0.5, 0.7))
        self.assertEqual(pl1.isvalid, True)
        self.assertEqual(pl1.iscontinuous, True)
        self.assertEqual(pl1.isconvex, False)
        self.assertEqual(pl1.iscw, True)
        self.assertEqual(pl1.isccw, False)
        pi = IntexPoint((0.7, 0.7))
        pe = IntexPoint((0.6, 0.4))
        # pl1.plot(also_plot=[pi, pe])
        pl1.add_internal_point(p=pi)
        self.assertEqual(pi in pl1.internal, True)
        self.assertEqual(pl1.is_internal_point(p=pi), True)
        self.assertEqual(pl1.is_external_point(p=pi), False)
        pl1.add_external_point(p=pe)
        self.assertEqual(pi in pl1.external, False)
        self.assertEqual(all([pl1.is_internal_point(p=x) for x in pl1.internal]), True)
        self.assertEqual(all([pl1.is_external_point(p=x) for x in pl1.external]), True)

    # -----------------------
    # - TRIANGLE TESTS -
    # -----------------------

class TestTriangulation(unittest.TestCase):

    def setUp(self):
        self.p1 = Point((0.5, 0.5))
        self.p2 = Point((0., 1.5))
        self.p3 = Point((-0.5, 0.5))
        self.p4 = Point((-0.5, 0.5))
        self.ls1 = LineSegment(self.p1, self.p2)
        self.ls2 = LineSegment(self.p2, self.p3)
        self.ls3 = LineSegment(self.p3, self.p1)
        self.ls4 = LineSegment(self.p4, self.p1)

    def test_triangle_creation(self):
        tri = Triangle([self.ls1, self.ls2, self.ls3, self.ls4])  # created false to test purging
        self.assertEqual(len(tri.segments), 3)
        self.assertEqual(tri.purged, True)
        self.assertEqual(tri.isvalid, True)

        # this is not correct, calling __init__ should raise an exception
        self.assertRaises(Exception, Triangle.__init__, (self.ls1, self.ls2))

        # non-valid triangle: points are is_collinear, area zero
        p1 = Point((0.5, 0.5))
        p2 = Point((0., 0.5))
        p3 = Point((-0.5, 0.5))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p1)
        self.assertRaises(Exception, Triangle.__init__, (ls1, ls2, ls3))

    def test_triangle_area(self):
        tri = Triangle([self.ls1, self.ls2, self.ls3])
        self.assertEqual(tri.area, 0.5)
        self.assertDictEqual(tri.detailled_area, {'delaunay': 0.5, 'ear clipping': 0.5, 'triangle': 0.5})

    def test_triangle_angles_and_quality(self):
        p1 = Point((0.5, 0.5))
        p2 = Point((0., 1.0))
        p3 = Point((-0.5, 0.5))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p1)
        tri = Triangle([ls1, ls2, ls3])
        self.assertAlmostEqual(tri.angles[0], 45, delta=EPS)
        self.assertAlmostEqual(tri.angles[1], 45, delta=EPS)
        self.assertAlmostEqual(tri.angles[2], 90, delta=EPS)
        self.assertAlmostEqual(tri.quality, 60, delta=EPS)

    def test_congruency(self):
        tri1 = Triangle([self.ls1, self.ls2, self.ls3])
        tri2 = tri1.closed_copy()
        tri2.move(x0=0, y0=0)
        self.assertTrue(tri1.congruent(tri2))
        self.assertTrue(tri2.congruent(tri1))
        tri2.move(x0=1, y0=0)
        self.assertFalse(tri1.congruent(tri2))
        self.assertFalse(tri2.congruent(tri1))

    def test_enclosed(self):
        tri1 = Triangle([self.ls1, self.ls2, self.ls3])
        p1 = Point((0.25, 0.55))
        p2 = Point((0., 0.9))
        p3 = Point((-0.25, 0.55))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p1)
        tri2 = Triangle([ls1, ls2, ls3])
        self.assertFalse(tri1.congruent(tri2))
        self.assertTrue(tri1.enclosed(tri2))
        self.assertFalse(tri2.enclosed(tri1))

    def test_overlap(self):
        tri1 = Triangle([self.ls1, self.ls2, self.ls3])
        tri2 = tri1.closed_copy()
        # moved up a little, there are intersection_point lines but no internal points
        tri2.move(x0=0, y0=0.5)
        self.assertTrue(tri1.overlap(tri2))
        tri3 = tri1.closed_copy()
        # rotated about endpoint. intersectiong lines, conincident points but no internal points
        tri3.rotate_about_xy(x0=0.5, y0=0.5, phi=-30)
        self.assertTrue(tri1.overlap(tri3))
        # moved after rotation. general case.
        tri3.move(x0=-0.1, y0=0.1)
        self.assertTrue(tri1.overlap(tri3))
        # moved to be disjunct
        tri3.move(x0=10, y0=10)
        self.assertFalse(tri1.overlap(tri3))
        self.assertTrue(tri1.disjunct(tri3))

        # rotated about endpoint.
        # no intersecting lines, no internal points, colocated and touching points, overlapping lines
        p1 = Point((0.5, 0.))
        p2 = Point((0., math.sqrt(3)/2.))
        p3 = Point((-0.5, 0.))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p1)
        tri = Triangle([ls1, ls2, ls3])
        tri4 = tri.closed_copy()
        tri4.rotate_about_xy(x0=0.5, y0=0., phi=-(tri.angles[0]))
        self.assertFalse(tri.overlap(tri4))
        self.assertFalse(tri4.overlap(tri))
        self.assertTrue(tri.touch(tri4))
        self.assertFalse(tri.disjunct(tri4))

        tri4 = tri.closed_copy()
        tri4.rotate_about_xy(x0=0.5, y0=0., phi=-(tri.angles[0]+0.01))
        self.assertFalse(tri.overlap(tri4))
        self.assertFalse(tri4.overlap(tri))
        self.assertFalse(tri.touch(tri4))
        self.assertTrue(tri.disjunct(tri4))

        tri4 = tri.closed_copy()
        tri4.rotate_about_xy(x0=0.5, y0=0., phi=-(tri.angles[0]-0.01))
        self.assertTrue(tri.overlap(tri4))
        self.assertTrue(tri4.overlap(tri))
        self.assertFalse(tri.touch(tri4))

        tri4.move(x0=0.1, y0=0)
        self.assertFalse(tri.overlap(tri4))
        self.assertFalse(tri4.overlap(tri))
        self.assertFalse(tri.touch(tri4))
        self.assertTrue(tri.disjunct(tri4))

    # -----------------------
    # - TRIANGULATION TESTS -
    # -----------------------

    def test_triangle_area1(self):
        p1 = Point((0.5, 0.5))
        p2 = Point((0., 1.5))
        p3 = Point((-0.5, 0.5))
        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p1)
        tri = Triangle([ls1, ls2, ls3])
        tri.rotate_about_origin(phi=30)
        for v in tri.detailled_area.values():
            self.assertAlmostEqual(v, 0.5, delta=0.000001)
        self.assertAlmostEqual(tri.area, 0.5, delta=0.000001)

    def test_triangulation_square(self):
        """ concave square """
        self.pl1 = square(a=1, x0=0.5, y0=0.5)
        self.pl1._segments[0].i._set_coords(coords=(0.5, 0.7))  # one of the points moved
        self.pl1.rotate_about_origin(phi=30)
        for v in self.pl1.detailled_area.values():
            self.assertAlmostEqual(v, 0.9, delta=0.000001)

    def test_triangulation_star4(self):
        """ four-leg star """
        ls1 = LineSegment(Point((0.5, 0.5)), Point((0., 1.5)))
        ls2 = LineSegment(Point((0., 1.5)), Point((-0.5, 0.5)))
        ls3 = LineSegment(Point((-0.5, 0.5)), Point((-1.5, 0.)))
        ls4 = LineSegment(Point((-1.5, 0.)), Point((-0.5, -0.5)))
        ls5 = LineSegment(Point((-0.5, -0.5)), Point((0., -1.5)))
        ls6 = LineSegment(Point((0., -1.5)), Point((0.5, -0.5)))
        ls7 = LineSegment(Point((0.5, -0.5)), Point((1.5, 0.)))
        ls8 = LineSegment(Point((1.5, 0.)), Point((0.5, 0.5)))
        self.pl1 = PolygonLine([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8])
        for v in self.pl1.detailled_area.values():
            self.assertAlmostEqual(v, 3.0, delta=0.000001)

    def test_triangulation_rectangle_collinear(self):
        """ rectangle with lots of is_collinear points, rotated """
        p1 = Point((0., 0.))
        p2 = Point((0.1, 0.))
        p3 = Point((0.2, 0.))
        p4 = Point((0.3, 0.))
        p5 = Point((0.4, 0.))

        p6 = Point((0.4, 0.1))
        p7 = Point((0.4, 0.2))
        p8 = Point((0.4, 0.3))
        p9 = Point((0.5, 0.4))
        p10 = Point((0.5, 0.6))

        p11 = Point((0.3, 0.6))
        p12 = Point((0., 0.6))

        p13 = Point((0., 0.55))
        p14 = Point((0., 0.45))
        p15 = Point((0., 0.35))
        p16 = Point((0., 0.25))
        p17 = Point((0., 0.15))
        p18 = Point((0., 0.05))

        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p4)
        ls4 = LineSegment(p4, p5)
        ls5 = LineSegment(p5, p6)
        ls6 = LineSegment(p6, p7)
        ls7 = LineSegment(p7, p8)
        ls8 = LineSegment(p8, p9)
        ls9 = LineSegment(p9, p10)
        ls10 = LineSegment(p10, p11)
        ls11 = LineSegment(p11, p12)
        ls12 = LineSegment(p12, p13)
        ls13 = LineSegment(p13, p14)
        ls14 = LineSegment(p14, p15)
        ls15 = LineSegment(p15, p16)
        ls16 = LineSegment(p16, p17)
        ls17 = LineSegment(p17, p18)
        ls18 = LineSegment(p18, p1)
        from collections import deque
        segs = deque([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8, ls9, ls10, ls11, ls12, ls13, ls14, ls15, ls16, ls17, ls18])
        for i in range(len(segs)):
            segs.rotate(1)
            pl = PolygonLine(list(segs))
            for v in pl.detailled_area.values():
                self.assertAlmostEqual(v, 0.265, delta=0.000001)

    def test_triangulation_U(self):
        """ U-shape with lots of is_collinear points, rotated """
        p1 = Point((0., 0.))
        p2 = Point((1, 0.))
        p3 = Point((2, 0.))
        p4 = Point((3, 0.))
        p5 = Point((4, 0.))

        p6 = Point((4, 1))
        p7 = Point((4, 2))
        p8 = Point((4, 3))
        p9 = Point((4, 4))

        p10 = Point((3, 4))

        p11 = Point((3, 3))
        p12 = Point((3, 2))
        p13 = Point((3, 1))

        p14 = Point((2, 1))
        p15 = Point((1, 1))

        p16 = Point((1, 2))
        p17 = Point((1, 3))
        p18 = Point((1, 4))

        p19 = Point((0, 4))

        p20 = Point((0, 3))
        p21 = Point((0, 2))
        p22 = Point((0, 1))

        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p4)
        ls4 = LineSegment(p4, p5)
        ls5 = LineSegment(p5, p6)
        ls6 = LineSegment(p6, p7)
        ls7 = LineSegment(p7, p8)
        ls8 = LineSegment(p8, p9)
        ls9 = LineSegment(p9, p10)
        ls10 = LineSegment(p10, p11)
        ls11 = LineSegment(p11, p12)
        ls12 = LineSegment(p12, p13)
        ls13 = LineSegment(p13, p14)
        ls14 = LineSegment(p14, p15)
        ls15 = LineSegment(p15, p16)
        ls16 = LineSegment(p16, p17)
        ls17 = LineSegment(p17, p18)
        ls18 = LineSegment(p18, p19)
        ls19 = LineSegment(p19, p20)
        ls20 = LineSegment(p20, p21)
        ls21 = LineSegment(p21, p22)
        ls22 = LineSegment(p22, p1)

        from collections import deque
        segs = deque([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8, ls9, ls10, ls11, ls12, ls13, ls14, ls15, ls16, ls17, ls18, ls19, ls20, ls21, ls22])
        for i in range(len(segs)):
            segs.rotate(1)
            pl = PolygonLine(list(segs))
            # with removing is_collinear points first
            for v in remove_collinear_points(pl).detailled_area.values():
                self.assertAlmostEqual(v, 10., delta=0.000001)
            # without removing the is_collinear points
            for v in pl.detailled_area.values():
                self.assertAlmostEqual(v, 10., delta=0.000001)

    def test_triangulation_star6(self):
        """ six-leg star """
        pass

    def test_triangulation_simple_U1(self):
        """ U-shape, rotated """
        p1 = Point((0., -0.1))
        p2 = Point((4, 0.4))
        p3 = Point((4.3, 4.3))
        p4 = Point((3.1, 3.6))
        p5 = Point((3, 1))
        p6 = Point((1, 0.5))
        p7 = Point((0.9, 4.2))
        p8 = Point((-0.1, 5.1))

        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p4)
        ls4 = LineSegment(p4, p5)
        ls5 = LineSegment(p5, p6)
        ls6 = LineSegment(p6, p7)
        ls7 = LineSegment(p7, p8)
        ls8 = LineSegment(p8, p1)

        from collections import deque
        segs = deque([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8])
        for i in range(len(segs)):
            segs.rotate(1)
            pl = PolygonLine(list(segs))
            for v in pl.detailled_area.values():
                self.assertAlmostEqual(v, 9.8, delta=0.000001)

    def test_triangulation_simple_U2(self):
        """ irregular U-shape, rotated """
        p1 = Point((0., -0.1))
        p2 = Point((4, 0.4))
        p3 = Point((7.3, 4.3))
        p4 = Point((5.1, 3.6))
        p5 = Point((3, 1))
        p6 = Point((1, 0.5))
        p7 = Point((-7., 14.2))
        p8 = Point((-8.1, 5.1))

        ls1 = LineSegment(p1, p2)
        ls2 = LineSegment(p2, p3)
        ls3 = LineSegment(p3, p4)
        ls4 = LineSegment(p4, p5)
        ls5 = LineSegment(p5, p6)
        ls6 = LineSegment(p6, p7)
        ls7 = LineSegment(p7, p8)
        ls8 = LineSegment(p8, p1)

        from collections import deque
        segs = deque([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8])

        for i in range(len(segs)):
            segs.rotate(1)
            pl = PolygonLine(list(segs))
            for k, v in pl.detailled_area.items():
                self.assertAlmostEqual(v, 55.83, delta=0.000001)

    def test_triangulation_hole(self):
        """ square with a hole """
        pass
