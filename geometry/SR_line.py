# -*- coding: utf-8 -*-

"""
This is a module to deal with 2D lines, segments.
A segment - a LineSegment object - is defined by two endpoints (i, j, both Point objects) and the section between these.
The direction of the Segment is i -> j.
The endpoints are part of the SegmentLine. 

Relative position of SegmentLines
# A: they are is_parallel. self.intersection_point(other) returns None
# 1: No intersection_point point.
# 3: they overlap partially.
# 5: they overlap totally meaning: both endpoints are co-located.
# 2: they touch, that is, one of the endpoints is co-located. This also means there is an intersection_point Point

# B: They intersect. One intersection_point point, self.intersection_point(other) returns a Point instance.
# 1: The intersection_point point (IP) is off the segment.
# 2: The segments touch (on one of the endpoints or on the segment).
# 3: The intersection_point point is on the segment, not on endpoint

"""


from math import degrees, sqrt, atan2, acos
from geometry.SR_point import Point, EMPTY, distance, EPS
from geometry import _plotting_available, plt
from cached_property import cached_property
import itertools
import pprint as pp
import copy


class LineSegment(object):
    """
    LineSegment object defined by two endpoints: i, j
    Direction is from i to j.
    """
    def __init__(self, *args):
        if not args:
            args = (Point(), Point())
        self._endpoints = args
        self.plot_style = None

    def __str__(self):
        if self.isvalid:
            return self.__repr__()
        else:
            return 'Non-valid instance of class %s' % self.__class__.__name__

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.i, self.j)

    def __eq__(self, other):
        """ equality: endpoints and direction match """
        if not self.isvalid and not other.isvalid:
            return False

        return (
            type(other) == type(self) and
            self.i == other.i and self.j == other.j
        )

    def loose_equal(self, other):
        """ equality, regardless of direction """
        if not self.isvalid and not other.isvalid:
            return False

        return (
            type(other) == type(self) and set(self.ij) == set(other.ij) and self.norm == other.norm
        )

    def __ne__(self, other):
        """ not equal """
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.ij)

    def _set_endpoints(self, endpoints):
        """ set the endpoints. use this to overwrite, set endpoints """
        assert all(p.isvalid for p in endpoints)
        self._endpoints = endpoints

    @cached_property
    def isvalid(self):
        """ validity check for LineSegments """
        if not (isinstance(self.i, Point) or isinstance(self.j, Point)):
            return False
        if not (self.i.isvalid and self.j.isvalid):
            return False
        if self.isnull():
            return True

        return True

    @property
    def i(self):
        """ Point object at the beginiing of the line segment """
        return self._endpoints[0]

    @property
    def j(self):
        """ Point object at the end of the line segment """
        return self._endpoints[1]

    @property
    def ij(self):
        """ Both endpoints of the line segment """
        return self.i, self.j

    @property
    def delta_x(self):
        """ x projection of the difference of the endpoints """
        return self.j.x - self.i.x

    @property
    def delta_y(self):
        """ y projection of the difference of the endpoints """
        return self.j.y - self.i.y

    @property
    def norm(self):
        """ length of the line segment """
        return sqrt(self.delta_x ** 2 + self.delta_y ** 2)

    @property
    def as_vector(self):
        """ The segment as vector from the origin """
        return self.delta_x, self.delta_y

    @property
    def unitvector(self):
        """ unit vector in the direction of the line segment """
        if not self.isnull():
            return tuple([x / self.norm for x in self.as_vector])
        else:
            return 0, 0

    def crossproduct(self, other):
        """
        Cross product. self and other are transformed to have their starting points in the origin
        """
        assert other.isvalid

        a = self.as_vector
        b = other.as_vector
        return a[0] * b[1] - a[1] * b[0]

    def dotproduct(self, other):
        """
        Dot product. self and other are transformed to have their starting points in the origin.
        """
        assert other.isvalid
        a = self.as_vector
        b = other.as_vector
        return a[0] * b[0] + a[1] * b[1]

    def isnull(self):
        """ checks if line segment is null, that is endpoints are at the same location """
        if self.i == self.j:
            return True
        else:
            return False

    @property
    def direction(self):
        """
        direction of the line in degrees, measured fro the x axis CCW
        0 is east, 90 is north and so on
        """
        _ret = degrees(atan2(self.delta_y, self.delta_x))
        if _ret < 0:
            return 360 + _ret
        else:
            return _ret

    def is_internal_point(self, p=EMPTY):
        """
        Tells if a point is within the box (degenerate case: 1D-Box, aka line) that has self as one of the diagonals
        This doesn't tell if it is on the line (except for the degenerate case), for that see is_point_on_line
        """
        assert p.isvalid

        _i = self.i
        _j = self.j

        if min(_i.x, _j.x) - EPS < p.x <= max(_i.x, _j.x) + EPS\
                and min(_i.y, _j.y) - EPS <= p.y <= max(_i.y, _j.y) + EPS:
            return True
        else:
            return False

    def is_point_on_line(self, p=EMPTY):
        """
        Tells if a point is on the line of the line segment.
        This doesn't tell if it is within the endpoints for that see is_internal_point
        """
        assert p.isvalid

        if p == self.i or p == self.j:
            return True

        # print(abs(self.cross_product(LineSegment(self.i, p))))
        if abs(self.crossproduct_commonstart(LineSegment(self.i, p))) < EPS:

            return True
        else:
            return False

    def is_point_to_right(self, p=EMPTY):
        """
        Tells if the point is to the right (when looking in the direction of the line)
        If the point is left, or on the line, it returns False
        """
        assert p.isvalid

        # _val = self.crossproduct_commonstart(other=LineSegment(self.i, p))
        _val = self.crossproduct(other=LineSegment(self.i, p))
        if _val < 0:
            return True
        elif _val > 0:
            return False
        else:  # on the line - False
            return False

    def is_point_to_left(self, p=EMPTY):
        """
        Tells if the point is to the left (when looking in the direction of the line)
        If the point is left, or on the line, it returns False
        """
        assert p.isvalid

        if self.is_point_on_line(p=p):
            return False

        else:
            return not self.is_point_to_right(p=p)

    def opposite_direction(self, other):
        """
        Tells if the LineSegments have opposite direction
        This doesn't tell anything about being is_collinear, overlapping etc.
        """
        assert other.isvalid

        if (self.direction + 180) % 360 == other.direction:
            return True
        else:
            return False

    def is_parallel(self, other):
        """
        This tells if two LineSegments are is_parallel
        This doesn't tell if they have the same or opposite direction
        """
        assert other.isvalid

        if self.crossproduct(other) == 0:
            return True
        else:
            return False

    def is_intersecting(self, other):
        """
        This tells if two LineSegments are is_parallel
        This doesn't tell if they have the same or opposite direction
        """
        assert other.isvalid

        return not self.is_parallel(other)

    def line_point_distance(self, p):
        """
        This calculated the distance between the line of LineSegment and a point
        This doen't tell if the closest point on the line of the LineSegment is actually
        in/on the segment
        """
        assert p.isvalid

        return abs((self.j.y - self.i.y) * p.x - (self.j.x - self.i.x) * p.y +
                   self.j.x * self.i.y - self.j.y * self.i.x) / self.norm

    def continuous(self, other):
        """
        Tells if two LineSegments touch and the the common point is
        the endpoint of the one and the starting point of the other one.
        """
        assert other.isvalid

        if self.touch(other=other):
            if self.j == other.i or self.i == other.j:  # end-begin or begin-end
                return True

        return False  # either don'T touch or end-end or begin-begin

    def touch(self, other):
        """ Tells if two LineSegments have only one of their endpoints common """
        assert other.isvalid

        # 1. possibility: the nodes are either equal or identical
        if len(set(self.ij).intersection(set(other.ij))) == 1:
            _ret = True
        else:
            _ret = False

        # 2. possibility: the nodes are pretty close to each other
        for A in self.ij:
            for B in other.ij:
                if distance(A, B) < EPS:
                    _ret = True

        return _ret

    def overlap(self, other):
        """
        Tells if two parallel lines overlap, that is, they are on the same line
        Returns a tuple. [0] is boolean, True if they overlap
        [1] is None if no overlap, a segment representing the overlap part otherwise.
        The direction of the segment is irrelevant.
        If they touch, they don't overlap.
        """
        assert other.isvalid

        # if not parallel we stop right away
        if not self.is_parallel(other):
            # print('not parallel')
            return False, None

        # are they collinear?
        if not self.is_collinear(other):
            # print('not collinear')
            return False, None

        # no overlap at all: no internal points, no touching - disjunct lines
        if not any([self.is_internal_point(x) for x in other.ij]) \
                and not any([other.is_internal_point(x) for x in self.ij]) \
                and not self.touch(other):
            # print('nothing')
            return False, None

        # total overlap
        if self.loose_equal(other):
            _ret = copy.deepcopy(self)
            # print('total overlap')
            return True, _ret

        # one endpoint touches, one not. This other is either inside or outside
        if self.touch(other):
            cp = set(self.ij).intersection(set(other.ij))  # common point
            if len(cp) != 1:
                raise Exception('should be one as they touch and are not loose equal. check the evaluation of touch')
            cp = cp.pop()
            # the other points
            otherpoints = [x for x in self.ij + other.ij if x != cp]
            # two linesegments are created. since these are collinear, if their direction is the same: there is overlap
            l1 = LineSegment(cp, otherpoints[0])
            l2 = LineSegment(cp, otherpoints[1])
            if l1.direction == l2.direction:  # the two segments are in the same direction and cp is common point
                # print('overlapping segment')
                return True, sorted([l1, l2], key=lambda x: x.norm)[0]
            else:
                # print('touch in one point')
                return True, cp  # the segments are in the opposite direction

        # both endpoints of one of the segments is internal for the other
        a, b = other, self
        for i in range(2):
            a, b = b, a
            if all([a.is_internal_point(x) for x in b.ij]):
                _ret = copy.copy(b)
                # print('one segment fully in another')
                return True, _ret

        # one endpoint of one of the segments is internal for the other segment
        a, b = other, self
        _inpoints = []
        for i in range(2):
            a, b = b, a
            _inpoints += [x for x in b.ij if a.is_internal_point(x)]
        # print('partially overlapping segments')
        return True, LineSegment(*sorted(_inpoints, key=lambda x: distance(Point((0, 0)), x)))

    def intersection_point(self, other):
        """
        This gives the coordinates of the point that lies both on the lines of LineSegments self and other
        This doesn't tell if the closest point is in/on the LineSegment.
        Returns None if the lines are is_parallel.
        Returns a Point instance if there is an intersection_point point.
        """
        assert other.isvalid

        if self.is_parallel(other):
            return None

        else:
            _x1 = self.i.x
            _y1 = self.i.y
            _x2 = self.j.x
            _y2 = self.j.y

            _x3 = other.i.x
            _y3 = other.i.y
            _x4 = other.j.x
            _y4 = other.j.y

            _den = float((_x1 - _x2) * (_y3 - _y4) - (_y1 - _y2) * (_x3 - _x4))
            _num_x = (_x1 * _y2 - _y1 * _x2) * (_x3 - _x4) - (_x1 - _x2) * (_x3 * _y4 - _y3 * _x4)
            _num_y = (_x1 * _y2 - _y1 * _x2) * (_y3 - _y4) - (_y1 - _y2) * (_x3 * _y4 - _y3 * _x4)
            return Point(((_num_x / _den), (_num_y / _den)))

    def common_part(self, other):
        """
        tells if there is a common part of two segments and what it is.
        intersecting lines:
        - intersection_point on both segments, inkc. Endpoints: Point
        - intersection_point off both segments: None
        
        parallel lines:
        - non is_collinear segments: None
        - is_collinear segments:
            - partial overlap: the overlapping part as SegmentLine
            - full overlap (also: two common endpoints): SegmentLine
            - 1 common endpoint (touching lines): Point
            - no common part: None
        """

        if self.is_parallel(other):
            # print('')
            # print('parallel')
            if self.is_collinear(other):
                # print('collinear')
                # full overlap, endpoints are equal, regardless of direction
                if self.loose_equal(other):
                    # print('loose equal')
                    _ret = copy.deepcopy(self)
                    return _ret

                else:  # touch etc. - all other cases
                    _ol = self.overlap(other)  # [0]: T/F, [1] the overlap as Segment
                    return _ol[1]

            else:
                return None

        else:
            ip = self.intersection_point(other)  # intersection point
            if self.intersection_on_segment(other) or self.touch(other):
                return ip
            else:
                return None

    def common_part_type(self, other):
        cp = self.common_part(other)
        if cp:
            return cp.__class__.__name__

        else:
            return None

    def is_collinear(self, other):
        """ tells if two lines are is_collinear """
        return all([self.is_point_on_line(x) for x in other.ij])

    def intersection_on_endpoint(self, other):
        """ This tells if there is an intersection_point and if the intersection_point point is one of the endpoints of 
        any of the endpoints of any of the segments.
        the point is not necessarily on _both_ segments
        """
        assert other.isvalid

        # calculating the intersection_point point of the lines
        _inters = self.intersection_point(other)

        if _inters is None:  # no intersection_point point - the lines are parallel
            return False

        # the intersection_point is a point where two endpoints are located
        if self.touch(other):
            return True

        # check if the intersection_point point is any of the endpoints
        else:
            if any([_inters == x for x in itertools.chain.from_iterable([self.ij, other.ij])]):
                return True
            else:
                return False

    def intersection_on_segment(self, other):
        """ This tells if there is an intersection_point and if it is on the segment, that is between the endpoints """
        assert other.isvalid

        if self.touch(other):
            return False

        _inters = self.intersection_point(other)

        if _inters is None:
            return False

        else:
            if self.is_internal_point(p=_inters) and other.is_internal_point(p=_inters):
                return True
            else:
                return False

    def mirror_vertical(self, x0=0):
        """ Mirrors the segment on a vertical line located at x0 """
        for _p in self._endpoints:
            _p._set_coords(coords=(_p.x - 2 * (_p.x - x0), _p.y))

    def move(self, x0=None, y0=None):
        """ Moves the line by x0, y0 units """
        for _p in self._endpoints:
            _p._set_coords(coords=(_p.x + x0, _p.y + y0))

    def reverse(self):
        """ reverses the segment """
        self._endpoints = self._endpoints[::-1]

    @property
    def line_reversed(self):
        """ Returns a LineSegment object with reserved endpoints """
        return LineSegment(self.j, self.i)

    # def angle_at_point(self):
    #     """ The angle between two segments that connect to the same point """
    #     A = acos(dot(v1, v2) / (v1.length() * v2.length()))

    def angle_between_lines(self, other):
        """ The angle between two LineSegment objects """
        assert other.isvalid

        # return atan2(norm(self.cross_product(other)), self.dot_product(other))
        _ret = degrees(acos(self.dotproduct(other) / (self.norm * other.norm)))
        return _ret

    def angle_from_start_point(self, p=EMPTY):
        """ The angle between self and the LineSegment defined by self.i and p """
        assert p.isvalid
        assert p != self.i

        return self.angle_between_lines(LineSegment(self.i, p))

    def crossproduct_commonstart(self, other):
        """
        https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python
        2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
        Returns a positive value, if OAB makes a counter-clockwise turn, (that is, right)
        negative for clockwise turn, (that is, left)
        and zero if the point_set are is_collinear.
        """

        assert other.isvalid
        assert self.i == other.i

        return (self.j.x - self.i.x) * (other.j.y - self.i.y) - (self.j.y - self.i.y) * (other.j.x - self.i.x)

    @property
    def bounding_box(self):
        return (min([z.x for z in self.ij]), max([z.x for z in self.ij])), (min([z.y for z in self.ij]), max([z.y for z in self.ij]))

    @property
    def midpoint(self):
        return (2 * self.i.x + self.j.x) / 3, (2 * self.i.y + self.j.y) / 3

    def plot(self, show=False, style=None, annotate=True):

        if _plotting_available:
            if style is None:
                self.plot_style = 'b-'
            else:
                self.plot_style = style

            plt.plot([p.x for p in self.ij], [p.y for p in self.ij], self.plot_style)
            mp = self.midpoint
            ax = plt.gca()
            ax.arrow(mp[0], mp[1], self.unitvector[0] * (self.norm / 3.), self.unitvector[1] * (self.norm / 3.),
                     head_width=self.norm / 30., head_length=self.norm / 15., fc='k', ec='k')
            if show:
                plt.axis('tight')
                plt.axis('equal')
                plt.show()


def point_on_line_at_ratio(line, ratio):
    """
    this provides a point that lies on the line at a given ratio
    the ratio is understood as the ratio to the distance from the endpoint i to the full length of the line
    0 - Point i
    1 - Point j
    """
    assert 0 <= ratio <= 1

    if ratio == 0:
        return line.i

    elif ratio == 1:
        return line.j

    else:
        linelength = line.norm
        _x = (1 - ratio) * line.i.x + ratio * line.j.x
        _y = (1 - ratio) * line.i.y + ratio * line.j.y
        return Point((_x, _y))


def divide_line(line, numseg):
    """ this divides a line in numseg segments of equal length """
    numseg = float(numseg)
    _pts = [line.i]
    for i in range(1, int(numseg)):
        _x = ((numseg - i) / numseg) * line.i.x + (i / numseg) * line.j.x
        _y = ((numseg - i) / numseg) * line.i.y + (i / numseg) * line.j.y
        _pts.append(Point((_x, _y)))
    _pts.append(line.j)

    _ret = []
    for p1, p2 in zip(_pts, _pts[1:]):
        _ret.append(LineSegment(p1, p2))
    return _ret


def add_point_on_line(line, p):
    """ this adds p to the line - divides the line - if the point is internal. returned are two linesegments """
    if line.is_internal_point(p=p):
        return LineSegment(line.i, p), LineSegment(p, line.j)
    else:
        return False


def run():
    pass


if __name__ == '__main__':
    run()
