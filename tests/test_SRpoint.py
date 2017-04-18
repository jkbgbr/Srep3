# -*- coding: utf-8 -*-


from geometry.SR_point import Point, IntexPoint, EMPTY
from geometry.SR_line import LineSegment
import unittest

# http://www.drdobbs.com/testing/unit-testing-with-python/240165163


class TestPoints(unittest.TestCase):
    """

    """

    def test_input_validity(self):
        p = Point()
        self.assertEqual(p.isvalid, False)
        p = Point((0, 1))
        self.assertEqual(p.isvalid, True)
        p = Point([(0, 1), (0, 1)])
        self.assertEqual(p.isvalid, False)
        p = Point([(0, (0, 1))])
        self.assertEqual(p.isvalid, False)
        p = Point([('s', (0, 1))])
        self.assertEqual(p.isvalid, False)
        p = Point([('s', 'mi')])
        self.assertEqual(p.isvalid, False)
        p = Point([(1, 's')])
        self.assertEqual(p.isvalid, False)
        p = Point([(1, 2, 3)])
        self.assertEqual(p.isvalid, False)
        self.assertEqual(Point(EMPTY).isvalid, False)

    def test_equality(self):
        self.assertEqual(Point((0, 0)), Point((0, 0)))
        self.assertNotEqual(Point((0, 1)), Point((0, 0)))
        # empty Points
        p1 = Point()
        p2 = Point((0, 0))
        p2._empty()
        self.assertNotEqual(p1, p2)
        # point set
        ps1 = [Point((0, 1)), Point((0, 0))]
        ps2 = [Point((0, 1)), Point((0, 0))]
        ps3 = [Point((0, 1)), Point((0, 1))]
        self.assertEqual(ps1, ps2)
        self.assertNotEqual(ps1, ps3)
        p1 = Point((0, 1))
        p2 = Point((0, 1))
        self.assertEqual(p1 == p2, True)
        p2 = Point((0, -1))
        self.assertEqual(p1 == p2, False)

    def test_magic_methods(self):
        p1 = Point((0, 1))
        p2 = Point((0, 1))
        self.assertEqual(p1+p2, Point((0, 2)))
        self.assertEqual(p1-(p2+p2), Point((0, -1)))

        p3 = Point([(0, 1), (0, 1)])
        p4 = Point([(-1, 21), (4, -71)])
        self.assertEqual(p3+p4, 'Non-valid instance of class Point')
        self.assertEqual(p3-p4, 'Non-valid instance of class Point')


    def test_coordinate_change(self):
        ps2 = [Point((0, 1)), Point((0, 0))]
        ps3 = [Point((0, 1)), Point((0, 1))]
        ps3[1]._set_coords(coords=(ps2[1].x, ps2[1].y))
        self.assertEqual(ps2, ps3)

        p1 = Point((0, 1))
        p2 = Point((0, 0))
        self.assertNotEqual(p1, p2)
        p2._set_coords(coords=(0, 1))
        self.assertEqual(p1, p2)
        p2._set_coords(coords=p1.xy)
        self.assertEqual(p1, p2)

        p1 = Point((1, 2))
        p1._empty()
        self.assertEqual(p1._isempty(), True)
        p1._set_coords(coords=(1, 3))
        self.assertEqual(p1._isempty(), False)

    def test_repr(self):
        p1 = Point((0, 1))
        self.assertEqual(p1.__str__(), 'Point((0.0, 1.0))')
        p1 = Point()
        self.assertEqual(p1.__str__(), 'Non-valid instance of class Point')

    def test_xy(self):
        p1 = Point((0, 1))
        self.assertEqual(p1.xy, (0, 1))
        p1 = Point((1, 2))
        self.assertEqual(p1.xy, (1, 2))
        p1 = Point((1, 2))
        p1._set_coords(coords=(1, 3))
        self.assertEqual(p1.xy, (1, 3))

    def test_point_distance(self):
        p1 = Point((0, 0))
        p2 = Point((3, 4))
        from geometry.SR_point import distance
        self.assertEqual(distance(p1, p2), 5)
        p1._set_coords(coords=(1, 3))
        self.assertEqual(p1._isempty(), False)

    def test_distance(self):
        p3 = Point((0, 5))
        p4 = Point((1, 5))
        self.assertEqual(p3.distance(p4), 1)

    def test_bounding_box(self):
        p3 = Point((3, 5))
        self.assertTupleEqual(p3.bounding_box(), ((3, 3), (5, 5)))

    def test_intex_point(self):
        p1 = Point((3, 5))
        p2 = Point((2, 9))
        s1 = LineSegment([p1, p2])
        p3 = IntexPoint((8, 1), segment=s1)


if __name__ == '__main__':
    unittest.main()
