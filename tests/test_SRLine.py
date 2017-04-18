# -*- coding: utf-8 -*-

from geometry.SR_line import LineSegment, divide_line, add_point_on_line, point_on_line_at_ratio
from geometry.SR_point import Point, EPS, PRECISION
import unittest
import math

SMALL_DISTANCE = 1./(10.**(PRECISION-5))


class TestLines(unittest.TestCase):
    """

    """

    def setUp(self):
        self.l0 = LineSegment(Point((0, 0)), Point((0, 0)))
        self.l1 = LineSegment(Point((0, 1)), Point((0, 2)))
        self.l2 = LineSegment(Point((0, 1)), Point((0, 2)))
        self.l3 = LineSegment(Point((0, 1)), Point((0, 3)))
        self.l4 = LineSegment(Point((0, 0)), Point((2, 2)))
        self.l5 = LineSegment(Point((0, 2)), Point((0, 1)))

    def test_repr(self):
        self.assertEqual(self.l1, eval(self.l1.__repr__()))
        self.assertEqual(self.l1.__str__(), 'LineSegment(Point((0.0, 1.0)), Point((0.0, 2.0)))')

    def test_if_null(self):
        # checking if LineSegment is null
        l = LineSegment()
        self.assertEqual(l.isnull(), False)
        self.assertEqual(self.l1.isnull(), False)

    def test_endpoint_Setting(self):
        import copy
        l1 = copy.deepcopy(self.l1)
        p1 = Point((0, 56))
        p2 = Point((48, 56))
        l1._set_endpoints([p1, p2])
        self.assertEqual(l1, LineSegment(p1, p2))
        self.assertEqual(l1.loose_equal(LineSegment(p2, p1)), True)
        self.assertEqual(l1.loose_equal(LineSegment(p1, p2)), True)

    def test_equality(self):
        # equality
        self.assertEqual(self.l1, self.l2)
        self.assertNotEqual(self.l1, self.l3)
        self.assertNotEqual(self.l1, self.l5)

        self.assertEqual(self.l1.ij, (self.l1.i, self.l1.j))

    def test_deltas(self):
        self.assertEqual(self.l1.delta_x, 0)
        self.assertEqual(self.l1.delta_y, 1)
        self.assertEqual(self.l4.delta_x, 2)
        self.assertEqual(self.l4.delta_y, 2)

    def test_norm(self):
        self.assertEqual(self.l1.norm, 1)
        self.assertEqual(self.l3.norm, 2)
        self.assertEqual(self.l4.norm, math.sqrt(8))

    def test_unitvector(self):
        self.assertTupleEqual(self.l1.unitvector, (0, 1))
        self.assertTupleEqual(self.l0.unitvector, (0, 0))
        self.assertTupleEqual(self.l4.unitvector, (1 / math.sqrt(2), 1 / math.sqrt(2)))

    def test_direction(self):
        l5 = LineSegment(Point((0, 0)), Point((1, 0)))
        l6 = LineSegment(Point((0, 0)), Point((1, 1)))
        l7 = LineSegment(Point((0, 0)), Point((0, 1)))
        l8 = LineSegment(Point((0, 0)), Point((-1, 1)))
        l9 = LineSegment(Point((0, 0)), Point((-1, 0)))
        l10 = LineSegment(Point((0, 0)), Point((-1, -1)))
        l11 = LineSegment(Point((0, 0)), Point((0, -1)))
        l12 = LineSegment(Point((0, 0)), Point((1, -1)))
        l13 = LineSegment(Point((0, 0)), Point((0, 0)))

        self.assertEqual(l5.direction, 0)
        self.assertEqual(l6.direction, 45)
        self.assertEqual(l7.direction, 90)
        self.assertEqual(l8.direction, 135)
        self.assertEqual(l9.direction, 180)
        self.assertEqual(l10.direction, 225)
        self.assertEqual(l11.direction, 270)
        self.assertEqual(l12.direction, 315)
        self.assertEqual(l13.direction, 0)

    def test_validity(self):
        # testing for validity
        l1 = LineSegment(Point((0, 0)), Point((0, 0)))
        self.assertEqual(l1.isvalid, True)
        l1 = LineSegment((0, 0), (0, 0))
        self.assertEqual(l1.isvalid, False)
        l1 = LineSegment(Point(), Point())
        self.assertEqual(l1.isvalid, False)
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.assertEqual(l1.isvalid, True)

    def test_is_internal_point(self):
        # checking if a point is internal
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.assertEqual(l1.is_internal_point(Point((0, 2))), False)
        self.assertEqual(l1.is_internal_point(Point((0, 0.5))), True)

    def test_point_on_line(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.assertEqual(l1.is_point_on_line(Point((1, 2))), False)
        self.assertEqual(l1.is_point_on_line(Point((0, 2))), True)
        self.assertEqual(l1.is_point_on_line(Point((0, 0.5))), True)
        self.assertEqual(l1.is_point_on_line(Point((0, 1))), True)
        self.assertEqual(l1.is_point_on_line(Point((0, -1))), True)

    def test_angle_from_start_point(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.assertEqual(l1.angle_from_start_point(Point((1, 0))), 90)
        # self.assertEqual(l1.angle_from_end_point(Point((1, 0))), 225)

    def test_line_point_distance(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        self.assertEqual(l1.line_point_distance(Point((1, 1))), 1)
        l1 = LineSegment(Point((0, 0)), Point((1, 1)))
        self.assertAlmostEqual(l1.line_point_distance(Point((0, 1))), math.sqrt(2) / 2, delta=EPS)

    def test_segments_parallel(self):
        l1 = LineSegment(Point((0, 0)), Point((1, 1)))
        l2 = LineSegment(Point((1, 2)), Point((2, 3)))
        self.assertEqual(l1.is_parallel(other=l2), True)
        l2 = LineSegment(Point((0, 0)), Point((1, 1)))
        self.assertEqual(l1.is_parallel(other=l2), True)
        l2 = LineSegment(Point((5, 5)), Point((8, 9)))
        self.assertEqual(l1.is_parallel(other=l2), False)
        self.assertEqual(l1.is_intersecting(other=l2), True)
        l2 = LineSegment(Point((5, 5)), Point((4, 4)))
        self.assertEqual(l1.is_parallel(other=l2), True)
        self.assertEqual(l1.is_intersecting(other=l2), False)

    def test_intersection_of_lines(self):
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((0, 2)), Point((2, 0)))
        self.assertEqual(l1.intersection_point(l2), Point((1, 1)))
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((1, 0)), Point((3, 2)))
        self.assertEqual(l1.intersection_point(l2), None)
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((3, 0)), Point((4, 1)))
        self.assertEqual(l1.intersection_point(l2), None)
        l1 = LineSegment(Point((0, 0)), Point((15, 4)))
        l2 = LineSegment(Point((0, 2)), Point((2, 0)))
        self.assertEqual(l1.intersection_point(l2), Point((1.5789473684210527, 0.42105263157894735)))

    def test_common_part(self):
        # parallel cases
        # not collinear
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((1, 0)), Point((3, 2)))
        self.assertIsNone(l1.common_part(l2))

        # loose equal A
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((2, 2)), Point((0, 0)))
        self.assertEqual(l1.common_part(l2), l1)

        # equal
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        l2 = LineSegment(Point((0, 0)), Point((2, 2)))
        self.assertEqual(l1.common_part(l2), l1)

        # 1a: both are internal, as l2 is inside l1
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((2, 2)), Point((3, 3)))
        self.assertEqual(l1.common_part(l2), l2)

        # 1b: both are internal, as l1 is inside l2
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((-2, -2)), Point((13, 13)))
        self.assertEqual(l1.common_part(l2), l1)

        # 2a: one endpoint touches and there is overlap
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((0, 0)), Point((3, 3)))
        self.assertEqual(l1.common_part(l2), l2)

        # 2b: one endpoint touches and there is NO overlap
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((0, 0)), Point((-3, -3)))
        self.assertEqual(l1.common_part(l2), Point((0, 0)))

        # 3a: both endpoints of l1 are inside l2
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((1, 1)), Point((2, 2)))
        self.assertEqual(l1.common_part(l2), l2)

        # 3b: both endpoints of l2 are inside l1
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((-1, -1)), Point((12, 12)))
        self.assertEqual(l1.common_part(l2), l1)

        # 4a: partial overlap, one point inside, one outside
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((-1, -1)), Point((3, 3)))
        self.assertEqual(l1.common_part(l2), LineSegment(Point((0, 0)), Point((3, 3))))

        # 4b: partial overlap, one point inside, one outside
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((4, 4)), Point((13, 13)))
        self.assertEqual(l1.common_part(l2), LineSegment(Point((4, 4)), Point((5, 5))))
        self.assertNotEqual(l1.common_part(l2), LineSegment(Point((5, 5)), Point((4, 4))))

        # non-parallel cases
        # intersection on the segments
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((0, 5)), Point((5, 0)))
        self.assertEqual(l1.common_part(l2), Point((2.5, 2.5)))

        # intersection on the endpoints
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((5, 5)), Point((8, 0)))
        self.assertEqual(l1.common_part(l2), Point((5, 5)))

        # intersection somewhere
        l1 = LineSegment(Point((0, 0)), Point((5, 5)))
        l2 = LineSegment(Point((6, 6)), Point((6, 10)))
        self.assertEqual(l1.common_part(l2), None)

    def test_common_part_type(self):
        # paralel but not collinear - no common part
        l1 = LineSegment(Point((0, 0)), Point((0, 5)))
        l2 = LineSegment(Point((1, 0)), Point((1, 6)))
        self.assertEqual(l1.common_part_type(l2), None)

        # total overlap - loose equal - common part is a LineSegment
        l1 = LineSegment(Point((0, 0)), Point((0, 5)))
        l2 = LineSegment(Point((0, 5)), Point((0, 0)))
        self.assertEqual(l1.common_part_type(l2), 'LineSegment')

        # intersection on line - common part is a Point
        l1 = LineSegment(Point((0, 0)), Point((0, 5)))
        l2 = LineSegment(Point((-1, 3)), Point((1, 4)))
        self.assertEqual(l1.common_part_type(l2), 'Point')



    def test_point_which_side(self):
        l1 = LineSegment(Point((0, 0)), Point((2, 2)))
        self.assertEqual(l1.is_point_to_right(p=Point((1, 1))), False)
        self.assertEqual(l1.is_point_to_right(p=Point((2, 0))), True)
        self.assertEqual(l1.is_point_to_right(p=Point((2.1, 1.9))), True)
        self.assertEqual(l1.is_point_to_right(p=Point((2.+SMALL_DISTANCE, 2-SMALL_DISTANCE))), True)
        self.assertEqual(l1.is_point_to_right(p=Point((2.-SMALL_DISTANCE, 2))), False)
        self.assertEqual(l1.is_point_to_right(p=Point((-2, 0))), False)

        self.assertEqual(l1.is_point_on_line(p=Point((1, 1))), True)
        self.assertEqual(l1.is_point_to_left(p=Point((1, 1))), False)
        self.assertEqual(l1.is_point_to_left(p=Point((2, 0))), False)
        self.assertEqual(l1.is_point_to_left(p=Point((2.1, 2))), False)
        self.assertEqual(l1.is_point_to_left(p=Point((1.5, 2))), True)
        self.assertEqual(l1.is_point_to_left(p=Point((-2, 0))), True)

        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        self.assertEqual(l1.is_point_to_right(p=Point((1, 1))), True)
        self.assertEqual(l1.is_point_to_right(p=Point((-2, 0))), False)
        self.assertEqual(l1.is_point_to_right(p=Point((SMALL_DISTANCE, 2))), True)
        self.assertEqual(l1.is_point_to_right(p=Point((-SMALL_DISTANCE, 2))), False)
        self.assertEqual(l1.is_point_to_right(p=Point((0, 3))), False)

        self.assertEqual(l1.is_point_to_left(p=Point((1, 1))), False)
        self.assertEqual(l1.is_point_to_left(p=Point((-2, 0))), True)
        self.assertEqual(l1.is_point_to_left(p=Point((SMALL_DISTANCE, 2))), False)
        # self.assertEqual(l1.is_point_to_left(p=Point((-EPS, 2))), True)
        self.assertEqual(l1.is_point_to_left(p=Point((0, 3))), False)

    def test_cross_product(self):
        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, 0)), Point((2, 0)))
        self.assertEqual(l1.crossproduct(l2), -4)
        self.assertEqual(l1.crossproduct_commonstart(l2), -4)
        self.assertEqual(l2.crossproduct_commonstart(l1), 4)
        self.assertEqual(l2.crossproduct_commonstart(l2), 0)

    def test_check_overlapping(self):
        # testing lines are overlapping
        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, 0)), Point((2, 0)))
        self.assertTupleEqual(l1.overlap(l2), (False, None))

        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, 3)), Point((0, 4)))
        self.assertTupleEqual(l1.overlap(l2), (False, None))

        l1 = LineSegment(Point((0, -99)), Point((0, 99)))
        l2 = LineSegment(Point((0, -4)), Point((0, -2)))
        # self.assertTupleEqual(l1.overlap(l2), (True, 2))
        self.assertTupleEqual(l1.overlap(l2), (True, LineSegment(Point((0.0, -4.0)), Point((0.0, -2.0)))))

        l1 = LineSegment(Point((0, 0)), Point((0, 99)))
        l2 = LineSegment(Point((0, -4)), Point((0, 101)))
        # self.assertTupleEqual(l1.overlap(l2), (True, 99))
        self.assertTupleEqual(l1.overlap(l2), (True, LineSegment(Point((0.0, 0.0)), Point((0.0, 99.0)))))

        l1 = LineSegment(Point((0, 0)), Point((0, 99)))
        l2 = LineSegment(Point((0, 80)), Point((0, 101)))
        # self.assertTupleEqual(l1.overlap(l2), (True, 19))
        self.assertTupleEqual(l1.overlap(l2), (True, LineSegment(Point((0.0, 80.0)), Point((0.0, 99.0)))))

        l1 = LineSegment(Point((0, 0)), Point((0, 99)))
        l2 = LineSegment(Point((0, 100)), Point((0, 101)))
        self.assertTupleEqual(l1.overlap(l2), (False, None))

        l1 = LineSegment(Point((0, 0)), Point((0, 9)))
        l2 = LineSegment(Point((0, 1)), Point((0, -2)))
        # self.assertTupleEqual(l1.overlap(l2), (True, 1))
        self.assertTupleEqual(l1.overlap(l2), (True, LineSegment(Point((0.0, 0.0)), Point((0.0, 1.0)))))

        # paralel but not collinear
        l1 = LineSegment(Point((0, 0)), Point((0, 5)))
        l2 = LineSegment(Point((1, 0)), Point((1, 6)))
        self.assertTupleEqual(l1.overlap(l2), (False, None))

        # total overlap - loose equal
        l1 = LineSegment(Point((0, 0)), Point((0, 5)))
        l2 = LineSegment(Point((0, 5)), Point((0, 0)))
        self.assertTupleEqual(l1.overlap(l2), (True, l1))

    def test_if_segments_of_opposite_direction(self):
        # testing opposite direction
        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, 0)), Point((0, -3)))
        self.assertEqual(l1.opposite_direction(l2), True)
        l2 = LineSegment(Point((0, 0)), Point((0, 5)))
        self.assertEqual(l1.opposite_direction(l2), False)
        l2 = LineSegment(Point((0, 0)), Point((1, 5)))
        self.assertEqual(l1.opposite_direction(l2), False)

    def test_angle_between_lines(self):
        # testing angle between lines
        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, -1)), Point((0, -3)))
        self.assertAlmostEqual(l1.angle_between_lines(l2), 180, delta=EPS)
        l1 = LineSegment(Point((0, 0)), Point((0, 1)))
        l2 = LineSegment(Point((0, 0)), Point((0, -1)))
        self.assertAlmostEqual(l1.angle_between_lines(l2), 180, delta=EPS)
        l2 = LineSegment(Point((0, 0)), Point((1, 1)))
        self.assertAlmostEqual(l1.angle_between_lines(l2), 45, delta=EPS)
        l2 = LineSegment(Point((0, 0)), Point((1, -1)))
        self.assertAlmostEqual(l1.angle_between_lines(l2), 135, delta=EPS)

    def test_continuity(self):
        # testing continuity
        l1 = LineSegment(Point((0, 0)), Point((0, 2)))
        l2 = LineSegment(Point((0, 2)), Point((0, -3)))
        self.assertEqual(l1.continuous(l2), True)
        l2 = LineSegment(Point((0, 0)), Point((0, 5)))
        self.assertEqual(l1.continuous(l2), False)
        l2 = LineSegment(Point((0, 1)), Point((1, 5)))
        self.assertEqual(l1.continuous(l2), False)

    def test_intersection(self):
        # testing intersection_point on segment
        # l1 and l2 overlap
        l1 = LineSegment(Point((0., 0.)), Point((0., 2.)))
        l2 = LineSegment(Point((0., 2.)), Point((0., -3.)))
        self.assertEqual(l1.intersection_on_segment(l2), False)

        # l4 ends on l3
        l3 = LineSegment(Point((-2.0, 2.0)), Point((4.0, 2.0)))
        l4 = LineSegment(Point((-0.5, 2.0)), Point((-0.5, 1.9)))
        self.assertEqual(l3.intersection_on_segment(l4), True)

        # l3.i == l4.j
        l3 = LineSegment(Point((-2.0, 2.0)), Point((4.0, 2.0)))
        l4 = LineSegment(Point((-2.0, 1.9)), Point((-2.0, 2.0)))
        self.assertEqual(l3.intersection_on_endpoint(l4), True)

        # l1 and l2 intersect
        l1 = LineSegment(Point((0., 0.)), Point((0., 2.)))  # same as previously
        l2 = LineSegment(Point((-1, 1)), Point((1, 2)))
        self.assertEqual(l1.intersection_on_segment(l2), True)

        # l1 and l2 intersect but not on the segment
        l2 = LineSegment(Point((-1, 10)), Point((1, 10)))
        self.assertEqual(l1.intersection_point(l2), Point((-0.0, 10.0)))
        self.assertEqual(l1.intersection_on_segment(l2), False)

        # l1 and l2 are collinear but do not overlap or touch
        l2 = LineSegment(Point((0, 3)), Point((0, 5)))
        self.assertEqual(l1.intersection_point(l2), None)
        self.assertEqual(l1.intersection_on_segment(l2), False)

        # l1 and l2 are parallel but not collinear
        l2 = LineSegment(Point((1, 3)), Point((1, 5)))
        self.assertEqual(l1.intersection_point(l2), None)
        self.assertEqual(l1.intersection_on_segment(l2), False)

        # l1.i == l2.i
        l2 = LineSegment(Point((0, 0)), Point((1, 5)))
        self.assertEqual(l1.intersection_point(l2), Point((0.0, 0.0)))
        self.assertEqual(l1.intersection_on_segment(l2), False)
        self.assertEqual(l1.intersection_on_endpoint(l2), True)
        self.assertEqual(l1.touch(l2), True)



    def test_reversal(self):
        # testing line reversal
        self.assertTupleEqual(self.l1.line_reversed.ij, LineSegment(self.l1.j, self.l1.i).ij)
        l1 = LineSegment(Point((0, 3)), Point((0, 5)))
        l2 = LineSegment(Point((0, 5)), Point((0, 3)))
        l1.reverse()  # reverses
        self.assertEqual(l1, l2)

    def test_move(self):
        # testing line move
        l1 = LineSegment(Point((0, 3)), Point((0, 5)))
        l2 = LineSegment(Point((1, 4)), Point((1, 6)))
        l1.move(x0=1, y0=1)  # reverses
        self.assertEqual(l1, l2)

    def test_midpoint(self):
        l1 = LineSegment(Point((0, 3)), Point((0, 5)))
        self.assertEqual(l1.midpoint, (0, 11/3.))

    def test_bounding_box(self):
        l1 = LineSegment(Point((-2, 3)), Point((0, 5)))
        self.assertEqual(l1.bounding_box, ((-2, 0), (3, 5)))

    def test_loose_equality(self):
        l1 = LineSegment(Point((0, 3)), Point((0, 5)))
        l2 = LineSegment(Point((0, 5)), Point((0, 3)))
        l3 = LineSegment(Point((0, 3)), Point((0, 5)))
        self.assertEqual(l1.loose_equal(l2), True)
        self.assertEqual(l1.loose_equal(l3), True)

    def test_line_division(self):
        # line division test
        self.line = LineSegment(Point((0, 0)), Point((4, 3)))
        _retval = tuple([LineSegment(Point((0.0, 0.0)), Point((1.33333333, 1.0))),
                         LineSegment(Point((1.33333333, 1.0)), Point((2.66666667, 2.0))),
                         LineSegment(Point((2.66666667, 2.0)), Point((4.0, 3.0)))])
        self.assertTupleEqual(_retval, tuple(divide_line(self.line, 3)), True)

    def test_line_division_by_ratio(self):
        _line = LineSegment(Point((0, 0)), Point((4, 3)))
        by_05 = Point((2, 1.5))
        self.assertEqual(point_on_line_at_ratio(line=_line, ratio=0.5), by_05)
        by_025 = Point((1, 0.75))
        self.assertEqual(point_on_line_at_ratio(line=_line, ratio=0.25), by_025)
        by_01 = Point((0.4, 0.3))
        self.assertEqual(point_on_line_at_ratio(line=_line, ratio=0.1), by_01)




    def test_insert_point(self):
        # Point addition test
        self.line = LineSegment(Point((0, 0)), Point((4, 3)))
        _retval = add_point_on_line(self.line, Point((2, 1.5)))
        self.assertTupleEqual((LineSegment(Point((0, 0)), Point((2, 1.5))), LineSegment(Point((2, 1.5)), Point((4, 3)))),
                              _retval, True)
