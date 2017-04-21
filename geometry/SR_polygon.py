# -*- coding: utf-8 -*-

import math
from geometry.SR_point import Point, IntexPoint, EMPTY, PRECISION, EPS, distance
from geometry.SR_line import LineSegment, divide_line
import copy
import itertools
from geometry import _plotting_available, plt, PURGE
import random
import pprint as pp
import collections
from cached_property import cached_property


TRIANGULATION_ALGORITHM = 'triangle'  # currently only triangle works
# TRIANGULATION_ALGORITHM = 'delaunay'
# possible values are: 'delaunay', 'ear clipping', 'triangle'


class PolygonLine(object):
    """
    PolygonLine object defined by at least two continuous LineSegment object
    """
    def __init__(self, args=None, polygongroup=None):
        if not args:
            args = (LineSegment(), LineSegment())
        self.polygongroup = polygongroup
        self._segments = list(args)
        self._internal_point = []
        self._indifferent_point = []
        self._external_point = []
        self.purged = False
        self._triangles = set()
        # self.delaunay = Delaunay2d(self)

    def __nonzero__(self):
        return True

    def __repr__(self):
        return 'PolygonLine(%r)' % self._segments

    def __eq__(self, other):
        """ Equal, if all segments are equal """
        return (
            type(other) == type(self) and
            len(self._segments) == len(other._segments) and
            all([s == o for s, o in zip(self._segments, other._segments)])
        )

    def loose_subset(self, other):
        """
        if other is a subset of self, based only on the endpoints. This means: some segments of other overlap
        some segments of self, but self has some segments that do not have a counterpart in other.
        True, if ALL endpoints of other is in the points of self or self and other are equal
        The direction of the segments or the order of these is not considered in this evaluation
        """
        selfpoints = set(itertools.chain.from_iterable([x.ij for x in self.segments]))
        otherpoints = set(itertools.chain.from_iterable([x.ij for x in other.segments]))
        if otherpoints.issubset(selfpoints):
            return True
        else:
            return False

    def loose_equal(self, other, direction_matters=False):
        """ equality, regardless of the order and regarding or regardless the direction of the segments """
        # no checking for validity, as order is not necessarily OK
        # if not self.isvalid and not other.isvalid:
        #     return False

        # to be loosely equal, they must have the same number of segments
        if len(self.segments) != len(other.segments):
            return False

        segment_loose_equality = False  # segment loose equality
        if direction_matters:  # direction matters
            if set(self.segments) == set(other.segments):
                segment_loose_equality = True

        else:  # whatever the direction
            # we must find for ALL segments in self a segment in other that is loosely equal
            segment_loose_equality = []
            for x in self.segments:  # each segment in self
                _eq = []
                for y in other.segments:  # each segment in other
                    # is there in other.segments, straight or reversed
                    _eq.append(any([z.loose_equal(x) for z in [y, y.line_reversed]]))
                segment_loose_equality.append(any(_eq))

        return all(segment_loose_equality)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((s.ij for s in self._segments))

    @property
    def triangles(self):
        return self._triangles

    def _delete_triangles(self):
        self._triangles = set()

    @property
    def segments(self):
        return self._segments

    def neighbour(self, p=EMPTY):
        _segs = self.find_segment_by_point(point=p, whichend='both')
        return [x.i for x in _segs if x.i != p] + [x.j for x in _segs if x.j != p]

    def find_segment_by_point(self, point, whichend='start'):
        """
        Finds a segment that belongs to a point.
        whichend tells, if this point should be a starting point or an endpoint
        returns a list
        """
        _seq = [s for s in self.segments if point in s.ij]  # at least one element
        # print(_seq)
        if not _seq:
            return []
        else:
            if whichend == 'start':
                _s = [x for x in _seq if point == x.i]
            elif whichend == 'end':
                _s = [x for x in _seq if point == x.j]
            elif whichend == 'both':
                _s = [x for x in _seq if point in x.ij]
            else:
                raise Exception('whichend either start, end, or both, not %s' % whichend)
            return _s

    @property
    def midpoint(self):
        """ the point at the average of the Points """
        _x0 = sum([z.x for z in self.point_set])
        _y0 = sum([z.y for z in self.point_set])
        return Point([x / len(self.point_set) for x in [_x0, _y0]])

    def disjunct(self, other):
        """
        two polygons are disjunct if none of their triangles (from triangulation) overlap or are congruent
        """
        _all_tri_pairs = itertools.combinations(itertools.chain.from_iterable([x.triangles for x in [self, other]]), 2)
        if any(x[0].overlap(x[1]) or x[0].congruent(x[1]) for x in _all_tri_pairs):
            return False
        else:
            return True

    def purge(self):
        """
        Reduces the polyline by searching for duplicates in the Point list and making sure these segments
        share a point (same object) rather than have coincident Points (equal objects).
        This function is explicitly called before all transformations.
        """

        if not PURGE:  # to disable "purging", set in __init__.py
            return False

        if self.purged:  # to speed things up, if purged, no need to do anything
            pass

        else:
            if len(self.point_set) == len(self.all_points):  # wasn't purged before but then it is not necessary
                self.purged = True

            else:
                # find duplicates and at later encounters set them equal with the first one
                seen = []  # list of duplicates
                for x in self.all_points:
                    _found = [y for y in seen if x == y]
                    if not _found:  # first encounter, store stuff
                        seen.append(x)
                    else:  # later encounters.
                        # As we may have branches, there are possibly more than one segments affected
                        assert len(_found) >= 1
                        _affected_segments = [s for s in self.segments if x in s.ij]
                        for _s in _affected_segments:
                            if _s.i == x:
                                _s._set_endpoints(endpoints=(_found[0], _s._endpoints[1]))
                            else:
                                _s._set_endpoints(endpoints=(_s._endpoints[0], _found[0]))
                            x = _found[0]
                self.purged = True

    def ordered_copy(self):
        import copy
        _ = copy.deepcopy(self)
        return _.order()

    def order(self):
        """
        Orders the segments of a PolygonLine to result a valid PolygonLine
        This works only with a cycle - the unordered segments must be a closed (but not continuous) segment
        If the segments provided do not constitue a valid polygon (e.g. not closed) then False is returned.
        If the segments provided can be used to create a valid polygon (e.g. closed, continuous) then True is returned.
        The created Polygon is _always_ the same, with the first segments endpoint i being the closest to the origin.
        note: should two such points exist the Polygon is obviously not valid
        """
        _akt = self._segments[0]
        _end = _akt.j
        _ret = []
        count = 0
        _stopsignal = False

        # # has the starting segment a next segment? lets check it first. This should speed up things a bit
        # _nxt_gen = (x for x in self._segments if (x.i == _end or x.j == _end) and x != _akt)
        # try:
        #     next(_nxt_gen)
        # except StopIteration:
        #     return False

        while not _stopsignal:
            count += 1
            if count > 2 * len(self._segments):
                _stopsignal = True

            # the next segment, that has its i endpoint at the j endpoint of _akt
            # if none found, we return a False
            _nxt_gen = (x for x in self._segments if x.i == _end and x != _akt)
            _flip = False
            try:
                _nxt = next(_nxt_gen)
            except StopIteration:
                _flip = True
                _nxt_gen = (x for x in self._segments if x.j == _end and x != _akt)
                try:
                    _nxt = next(_nxt_gen)
                except StopIteration:
                    return False

            # ok so there is a next segment
            if _nxt not in _ret:
                if _flip:
                    _nxt.reverse()
                _ret.append(_nxt)

            # the new actual segment and ist endpoint
            _end = _nxt.j
            _akt = _nxt

            if count > 500:
                raise Exception('WTF')

        # finished ordering the polygon
        # making sure the resulting polygon is always the same: sorting the segments so that the first is the one
        # that has its i endpoint the closest to the origin
        _first = sorted(_ret, key=lambda x: (distance(Point((0, 0)), x.i), x.i.x, x.i.y, x.norm))[0]
        _ret = collections.deque(_ret)  # making it a deque to be able to ratete it
        while _ret[0] != _first:
            _ret.rotate()  # rotate to the first

        self._segments = list(_ret)

        return True

    def append_point_at_end(self, p=EMPTY):
        """
        Adds a point to the end, that is, a new segment with the starting point previous endpoint
        This operation maintains the state of being purged or not, thus self.purged remains unchanged
        """
        assert p.isvalid  # a Point object

        self._segments.append(LineSegment(self.last_point, p))

    def append_point_at_beginning(self, p=EMPTY):
        """
        Adds a point to the beginning, that is, a new segment with the starting point provided here and 
        endpoint the former startpoint
        This operation maintains the state of being purged or not, thus self.purged remains unchanged
        """
        assert p.isvalid  # a Point object

        self._segments.insert(0, LineSegment(p, self.first_point))

    def add_internal_point(self, p=EMPTY, segment=None):
        """ add a point that is for sure on the inside of the polygon """
        self._internal_point.append(IntexPoint(args=p.xy, segment=segment))

    def add_indifferent_point(self, p=EMPTY, segment=None):
        """ like internal or external points, when it is not important whether in or out. Bei Flachboden. """
        self._indifferent_point.append(IntexPoint(args=p.xy, segment=segment))

    def add_external_point(self, p=EMPTY, segment=None):
        """ add a point that is for sure on the external of the polygon """
        self._external_point.append(IntexPoint(args=p.xy, segment=segment))

    @property
    def internal(self):
        return self._internal_point

    @property
    def indifferent(self):
        return self._indifferent_point

    @property
    def external(self):
        return self._external_point

    @property
    def all_points(self):
        """ All points in the polyline. Sorted list, first by the x, than by the y coordinate """
        return sorted(list(itertools.chain(*[s.ij for s in self.segments])), key=lambda z: (z.x, z.y))

    @property
    def point_set(self):
        """ 'Set' of points in the polyline. Sorted list, first by the x, than by the y coordinate """
        return sorted(list(set(self.all_points)), key=lambda z: (z.x, z.y))
        # return sorted(list(set(itertools.chain(*[s.ij for s in self.segments]))), key=lambda z: (z.x, z.y))

    @property
    def crossing_lines(self):
        """ gives a list of tuples of lines that cross each other """
        cases = itertools.permutations(self.segments, r=2)
        return [x for x in cases if x[0].intersection_on_segment(x[1])]

    @property
    def isvalid(self):
        if not all(isinstance(x, LineSegment) for x in self.segments):
            # print('not all segments LineSegment')
            return False
        if not all(x.isvalid for x in self.segments):
            # print('not all segments valid LineSegments')
            return False
        if not self.iscontinuous:
            # print('not continuous')
            return False
        if len(self.segments) < 1:
            # print('too short')
            return False
        if self.crossing_lines:
            # print('crossing lines')
            return False

        return True

    @property
    def discontinuity(self):
        _ret = []
        for x, y in zip(self.segments, self.segments[1:]):
            if not x.continuous(y):
                _ret.append(y)
        return _ret

    @property
    def iscontinuous(self):
        """ This checks if the polyline is continuous """
        # doesn't check for validity, as continuity is necessary for validity

        return all(x.continuous(y) for x, y in zip(self.segments, self.segments[1:]))

    def fix(self):
        """
        This fixes a polygonline.
        Given a closed but not continuous polygonline it starts somewhere and reverses the "next"
        segments until it becomes continuous.
        """
        raise NotImplementedError

    def next_segment(self, seg, isclosed=True):
        """
        Yields the segment _after_ if possible
        isclosed=True is strict: only closed, valid polygons
        isclosed=False: loose. Tries to find the segment, that seg.j = result.i.
        """

        if isclosed:
            assert self.isclosed

            _ret = [x for x in self.segments if x.i == seg.j]
            if len(_ret) != 0:
                return _ret[0]
            else:
                raise Exception('WTF')

        else:
            raise NotImplementedError

    @property
    def isclosed(self):
        """ This tells if the polygon is closed """
        # validity is not checked, see iscontinuous

        if self.iscontinuous and self.first_point == self.last_point:
            return True
        else:
            return False

    @property
    def first_point(self):
        """ First point of the polyline. If closed, the first Point of the first Segment defined """
        return self.first_segment.i

    @property
    def last_point(self):
        """
        Last point of the polyline. If closed, the last Point of the last Segment defined (and equals the first point)
        """
        return self.last_segment.j

    @property
    def first_segment(self):
        """ First segment of the polyline. If closed, the first Point of the first Segment defined """
        return self.segments[0]

    @property
    def last_segment(self):
        """
        Last segment of the polyline. If closed, the last Point of the last Segment defined (and equals the first point)
        """
        return self.segments[-1]

    @property
    def isopen(self):
        """ This tells if it is open """
        # validity is not checked

        return not self.isclosed

    def closed_copy(self):
        """ returns a copy of the polyline that is closed """
        if self.isclosed:
            return copy.deepcopy(self)
        else:
            _ret = copy.deepcopy(self)
            _ret.close()
            return _ret

    def close(self):
        """ adds a segment to close the plygon """
        if self.isopen:
            self._segments.append(LineSegment(self.last_point, self.first_point))
        else:
            pass

    @property
    def cps(self):
        """ Provides a generator of cross product values for further use """
        # no asserting

        _s = self.segments  # original list
        _sl = _s[1:]
        _sl.append(self.segments[0])  # list with the starting point shifted

        return (x.crossproduct(y) for x, y in zip(_s, _sl))

    @property
    def isconvex(self):
        """ Tells if the polygon is convex """
        assert self.isvalid

        # criterium is: sign of cross product of consequtive segments is the same
        if all([x >= 0 for x in self.cps]) or all([x <= 0 for x in self.cps]):  # Zero: 180 degree
            return True
        else:
            return False

    @property
    def iscw(self):
        """ Tells if point order is cw """
        if self.isconvex:  # pure positive or negative values
            if sum(self.cps) > 0:  # positive means ccw
                return False
            else:
                return True

        else:  # positive and negative values mixed
            _ps = len([x for x in self.cps if x >= 0])
            _ns = len([x for x in self.cps if x <= 0])
            if _ps > _ns:
                return False
            elif _ps < _ns:
                return True
            else:
                print(_ps)
                print(_ns)
                raise Exception("This shouldn't happen")

    @property
    def isccw(self):
        """ Tells if point order is ccw """
        return not self.iscw

    def is_internal_point(self, p=Point(), plotit=False):
        """ tells if given point is internal. True if internal """
        # assert self.isvalid
        assert p.isvalid

        ctrlpkt = 3

        # if the point is one of the Points of the polygon - False.
        if p in self.point_set:
            return False

        _bb = self.bounding_box()
        _x, _y = _bb[0][1], _bb[1][1]
        _finished = False
        _inters = []
        while not _finished:
            for i in range(ctrlpkt):
                # # old solution - a Point somewhere to the north east is taken
                P = Point((random.uniform(_x, 500 * _x), random.uniform(_y, 30 * _y)))
                ls = LineSegment(p, P)  # from p draw line in random direction

                # # new solution: a Line towards the center of gravity and then to infinity
                # _np = len(self.point_set)
                # _cg = Point([sum([p.x for p in self.point_set]) / _np, sum([p.y for p in self.point_set]) / _np])
                # ls = LineSegment(p, _cg)
                # ls.j._set_coords(coords=[ls.j.x + 10000 * ls.norm * ls.delta_x, ls.j.y + 10000 * ls.norm * ls.delta_y])

                if plotit:
                    ls.plot(show=False)
                    # this would plot the intersection_point point
                    # for x in self.segments:
                    #     if ls.intersection_on_segment(x):
                    #         ls.intersection_point(x).plot(show=False, style='ko')

                _inters.append(len([x for x in set(self.segments) if ls.intersection_on_segment(x)]))
            # check: all attempts should give the same result - instead of making sure the projected line is not on an
            # endpoint of one of the segments
            if not all([_inters[x] % 2 == _inters[0] % 2 for x in range(ctrlpkt)]):
                _finished = False
                _inters = []
            else:
                _finished = True

        # odd number of intersections means the point is INSIDE
        if _inters[0] % 2 == 1:
            return True
        else:
            return False

    def is_external_point(self, p=EMPTY):
        return not self.is_internal_point(p=p)

    def is_on_boundary(self, p=EMPTY):
        """ Tells if the point is on the Polygonline """
        if self.is_internal_point(p=p) or self.is_external_point(p=p):
            return False
        if p in self.point_set:
            return True
        if any([x.is_internal_point(p=p) and x.is_point_on_line(p=p) for x in self.segments]):
            return True
        return False

    def bounding_box(self):
        """ returns a tuple of tuples: (xmin, xmax), (ymin, ymax) """
        _xmin = 10e10
        _xmax = -10e10
        _ymin = 10e10
        _ymax = -10e10
        for _s in self.segments:
            _xmin = min(_xmin, min([ep.x for ep in _s.ij]))
            _xmax = max(_xmax, max([ep.x for ep in _s.ij]))
            _ymin = min(_ymin, min([ep.y for ep in _s.ij]))
            _ymax = max(_ymax, max([ep.y for ep in _s.ij]))

        return (_xmin, _xmax), (_ymin, _ymax)

    def rotate_about_origin(self, phi=0):
        """
        Rotates polygon by phi degrees about the origin of the global coordinate system.
        This may introduce some numerical errors, thus are the results rounded.
        """
        _cp = math.cos(math.radians(phi))  # sin(phi) shorthand
        _sp = math.sin(math.radians(phi))  # cos(phi) shorthand
        self.purge()
        for _p in self.point_set:
            _c = (_p.x * _cp - _p.y * _sp, _p.x * _sp + _p.y * _cp)  # "exact values"
            _c = tuple([round(x, int(PRECISION)) for x in _c])  # rounding to "PRRECISION" number of decimals
            _p._set_coords(coords=_c)

    def mirror_vertical(self, x0):
        """ Mirrors the segment on a vertical line located at x0 """
        self.purge()
        for _p in self.point_set:
            _p._set_coords(coords=(_p.x - 2 * (_p.x - x0), _p.y))

    def move(self, x0=None, y0=None):
        """ Moves the polygon by x0, y0 units """
        self.purge()
        # for s in self._segments:
        #     s.move(x0=x0, y0=y0)
        for _p in self.point_set:
            _p._set_coords(coords=(_p.x + x0, _p.y + y0))

    def rotate_about_xy(self, x0=None, y0=None, phi=None):
        """ Rotates polygon by phi degrees about given x, y point """
        self.move(x0=-x0, y0=-y0)  # move so that he origin is in x0, y0
        self.rotate_about_origin(phi=phi)  # rotate
        self.move(x0=x0, y0=y0)  # move back

    @property
    def is_triangle(self):
        if len(self.point_set) == 3 and len(self.segments) == 3 and self.isclosed and self.area > 0:
            return True
        else:
            return False

    @cached_property
    def area(self):
        """
        checks for the results of three different ways to calculate the area and yields one
        """
        if len(set(self.detailled_area.values())) == 1:
            return self.detailled_area[TRIANGULATION_ALGORITHM]

        else:
            pp.pprint(self.detailled_area)
            import time
            time.sleep(0.5)

            raise Exception('we have a problem here')

    @property
    def detailled_area(self):
        """ calculates the detailled_area of the Polygon using different algorithms. For this triangulation is used. """
        _ret = {}
        # for alg in ['delaunay', 'ear clipping']:
        # for alg in ['triangle']:
        # for alg in ['delaunay', 'ear clipping', 'triangle']:
        for alg in [TRIANGULATION_ALGORITHM]:
            self.triangulate(algorithm=alg, force=True)
            _ret[alg] = round(sum([x.area for x in self.triangles]), int(PRECISION))
        return _ret

    @property
    def collinear_points(self):
        """
        provides two sets of points based on the Points in segemnt.
        any number of points may line up to be is_collinear; sides are the points at the
        ends of such lines, mids are the midpoints. no distinction based on position whasoever.
        _mids are the midpoints, that is non-side points of the collinears
        _sides are the points that are at the end of lines containing the is_collinear points
        """

        # todo: give each element in sides a number telling how long the is_collinear line is
        _sidescount = {}
        _sides = set()  # points that my be on the side of a is_collinear situation
        _mids = set()  # Points that are midpoints
        for p in [x.i for x in self.segments]:
            nb = self.neighbour(p=p)
            _l = LineSegment(*nb)  # the line connecting the neighbours
            if _l.is_point_on_line(p=p):  # points are is_collinear
                _mids.add(p)  # add to mids - p is the middle point
                _sides.add(nb[0])  # add sides to sides
                _sides.add(nb[1])  # add sides to sides
                _sides = {x for x in _sides if x not in _mids}  # sides and mids are kept disjunct

                # todo: fix this stuff to be more robust
                # counting the importance of a given sides element
                for s in _sides:
                    if s in _sidescount.keys():
                        _sidescount[s] += 1
                    else:
                        _sidescount[s] = 0
                for m in _mids:
                    if m in _sidescount.keys():
                        _sidescount.pop(m)

        return _mids, _sides, _sidescount

    def triangulate(self, algorithm=TRIANGULATION_ALGORITHM, force=False):
        """
        wraps all triangluation algorithms. If force == True, a new triangulation will be calculated regardless
        if it was previously calculated or not.
        """
        if algorithm is None:
            raise Exception('No algorithm is provided')

        if algorithm == 'random':
            algorithm = random.choice(['delaunay', 'ear clipping', 'triangle'])

        if not self.triangles or force:  # if previously not triangulated or is forced to recalculate
            if algorithm is not None:
                if algorithm == 'delaunay':
                    self.triangulate_delaunay()
                elif algorithm == 'ear clipping':
                    self.ear_clipping()
                elif algorithm == 'triangle':
                    self.triangle_mesh()
                else:
                    raise Exception('Triangulation algorithm %s unknown' % algorithm)

    def triangulate_delaunay(self):
        """
        creates a Delaunay triangulation
        based on: http://code.activestate.com/recipes/579021-delaunay-triangulation/
        """

        self._delete_triangles()
        # if self.triangles:
        #     return True

        def validate(polygon, _triangle):
            # remove triangles that invalidate the constrains
            _ret = True
            for s in _triangle.segments:
                if any([s.intersection_on_segment(other=w) for w in polygon.segments]):
                    _ret = False
            if not polygon.is_internal_point(p=_triangle.midpoint):
                _ret = False
            if _triangle.area <= 0:
                _ret = False

            return _ret

        self.delaunay.triangulate()  # do the triangulation

        # defining the triangles for self and removing non-internal triangles
        for tri in self.delaunay.triangles:
            _tris = tri + [tri[0]]
            _pol = polygon_from_nodes([self.delaunay.points[x].xy for x in _tris])
            _candidate = Triangle(_pol.segments)

            if validate(self, _candidate):
                self._triangles.add(_candidate)

        _remained_points = [x for x in self.point_set if x not in self.delaunay.all_points]
        _remained_segments = set([x for x in self.segments if not any([x.loose_equal(y) for y in self.delaunay.all_segments])])

        # this is a try:
        # 1: for each point find the segments that contain it, and make a triangle from them. (polygon_from_nodes)
        # 2: add this triangle to polygon._triangles after validation
        # 3: repeat until no remainders (success!), OR remainders but no new triangles found/added (failure)

        # self.plot(show=False)
        # for tri in self.triangles:
        #     tri.plot(show=False)
        # tri.plot(show=True)

        while len(_remained_points) > 0:
            rp = _remained_points.pop()
            _segs = [x for x in _remained_segments if rp in x.ij]
            if len(_segs) == 2:
                _mi = list(itertools.chain.from_iterable(x.ij for x in _segs))
                _mi.append(rp)
                _mi = list(set(_mi))
                _mi = polygon_from_nodes(nodes=[w.xy for w in _mi])
                _mi.close()
                self._triangles.add(Triangle(_mi._segments))
            else:
                raise Exception('WTF')

    def ear_clipping(self):
        """
        creates a set of triangular polygons from the original polygon so that the triangles are not overlapping
        and only the points of the polygon are used.
        this is NOT a Delaunay triangulation and may fail.
        """

        self._delete_triangles()
        # if self.triangles:
        #     return True

        def _plot(_c_, _l_, _tri):
            # this plots the triangle
            _c_.plot(show=False)
            _l_.plot(show=False)
            _tri.midpoint.plot(show=True, style='k.')

        def chosen_point_ok(line, point):
            # checks if given point and line are is_collinear.
            if line.is_point_on_line(p=point):  # points are is_collinear
                return False
            else:
                return True

        def triangle_ok(polygon, triangle):
            # check if the midpoint of triangle is internal point to the polygon
            if not polygon.is_internal_point(p=triangle.midpoint):
                # print('nem belso haromszog')
                return False

            # check if any of the points of the triangle is external to the polygon
            if any([triangle.is_internal_point(p=x) for x in polygon.point_set if x not in triangle.point_set]):
                # print('belelognak')
                return False
            pairs = itertools.combinations(polygon.segments + triangle.segments, 2)
            for pair in pairs:
                if pair[0].intersection_on_segment(other=pair[1]):
                    # pair[0].plot(show=False)
                    # pair[0].plot(show=True)
                    return False

            # check if the area of the triangle is > 0
            if triangle.area <= 0:
                return False

            return True

        _c = copy.deepcopy(self)

        # _c.plot(show=True)

        # _c.order()
        count = 0
        cp_index = 0
        while len(_c.point_set) >= 3:
            count += 1

            if count > 5000:
                return False

            mids, sides, scount = _c.collinear_points
            sides = sorted(list(sides), key=lambda x: scount[x], reverse=True)

            try:
                cp = sides[cp_index]  # CurrentPoint
            except IndexError:
                cp = random.choice(tuple(_c.point_set))

            nb = _c.neighbour(p=cp)  # neightbours of cp

            # check if the chosen point and its two neighbours are is_collinear
            _l = LineSegment(*nb)  # the line connecting the neighbours
            if not chosen_point_ok(_l, cp):  # the chose point may not be is_collinear with the new segment
                cp_index += 1
                # print('point nem ok')
                continue

            # the new line defines a triange, this must have its midpoint inside the polygon
            _segs = _c.find_segment_by_point(point=cp, whichend='both') + [_l]
            tri = Triangle(_segs)  # this is the triangle defined by the neighbour points
            if not triangle_ok(polygon=_c, triangle=tri):
                cp_index += 1
                # print('triangle nem ok')
                continue

            # _plot(_c, _l, tri)
            self._triangles.add(tri)  # the triangle is added to the set

            # the original polygon is modified
            _c._segments = [x for x in _c._segments if x not in tri.segments]
            _c._segments += [_l]
            _c.purge()

        # for t in self.triangles:
        #     t.plot(show=False)
        #     t.midpoint.plot(show=False)
        # t.plot(show=True)

        return True

    def triangle_mesh(self):
        """
        this generates a constrined delaunay triangulation using the triangle lib and its python wrapper
        original lib: http://www.cs.cmu.edu/~quake/triangle.html
        wrapper: http://dzhelil.info/triangle/
        this requires numpy, matplotlib and triangle installed
        
        Perform triangulation on the input data `tri`. `tri` must be a dictionary
        that contains the following keys:
    
            * `vertices` - 2-D array that stores the xy position of each vertex
            * `segments` - optional 2-D array that stores segments. Segments are edges whose presence in the triangulation is enforced (although each segment may be subdivided into smaller edges). Each segment is specified by listing the indices of its two endpoints.
            * `holes` - optional 2-D array that stores holes. Holes are specified by identifying a point inside each hole. After the triangulation is formed, Triangle creates holes by eating triangles, spreading out from each hole point until its progress is blocked by PSLG segments; you must be careful to enclose each hole in segments, or your whole triangulation might be eaten away. If the two triangles abutting a segment are eaten, the segment itself is also eaten. Do not place a hole directly on a segment; if you do, Triangle will choose one side of the segment arbitrarily.
            * `regions` - optional 2-D array that stores region attributes and areas.
    
        The second (optional) arguments lists the options that should be passed to triangle.
    
            * `p` - Triangulates a Planar Straight Line Graph.
            * `r` - Refines a previously generated mesh.
            * `q` - Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q`.
            * `a` - Imposes a maximum triangle detailled_area constraint. A fixed detailled_area constraint (that applies to every triangle) may be specified after the `a`, or varying areas may be read from the input dictionary.
            * `c` - Encloses the convex hull with segments.
            * `D` - Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.
            * `X` - Suppresses exact arithmetic.
            * `S` - Specifies the maximum number of added Steiner points.
            * `i` - Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
            * `F` - Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
            * `l` - Uses only vertical cuts in the divide-and-conquer algorithm. By default, Triangle uses alternating vertical and horizontal cuts, which usually improve the speed except with vertex sets that are small or short and wide. This switch is primarily of theoretical interest.
            * `s` - Specifies that segments should be forced into the triangulation by recursively splitting them at their midpoints, rather than by generating a constrained Delaunay triangulation. Segment splitting is true to Ruppert's original algorithm, but can create needlessly small triangles. This switch is primarily of theoretical interest.
            * `C` - Check the consistency of the final mesh. Uses exact arithmetic for checking, even if the -X switch is used. Useful if you suspect Triangle is buggy.
        
        """
        import numpy as np
        import triangle

        # import matplotlib.pyplot as plt

        def triangle_ok(polygon, triangle):
            """ essentially the same as triangle_ok in the ear_clipping algorithm """
            # check if the midpoint of triangle is internal point to the polygon
            if not polygon.is_internal_point(p=triangle.midpoint):
                # print('midpoint is not internal!')
                return False
            # check if any of the points of the triangle is external to the polygon
            if any([triangle.is_internal_point(p=x) for x in polygon.point_set if x not in triangle.point_set]):
                # print('has external point!')
                return False

            # pairs = itertools.combinations(polygon.segments + triangle.segments, 2)
            # for pair in pairs:
            #     if not pair[0].intersection_on_endpoint(other=pair[1]):
            #         return False

            # check if the detailled_area of the triangle is > 0
            if triangle.area <= 0:
                print('A <= 0!')
                return False

            return True

        self._delete_triangles()

        # vertices are stored in a dict to have the connection between node numbering and node coordinates
        # here not only the coordinates are stored but thew Point instance to account for float inaccuracies when
        # evaluating equality
        # {0: Point((-1.0, 7.0)),
        #  1: Point((-1.0, 9.0))}
        _ap = self.point_set
        vertices = {k: v for k, v in zip(range(len(_ap)), _ap)}

        # segments
        # {Point((-1.0, 9.0)): 1,
        #  Point((-1.0, 7.0)): 0}
        _segs = {v: k for k, v in vertices.items()}
        # [(1, 0)]
        segments = [(_segs[s.i], _segs[s.j]) for s in self.segments]

        # data for triangle
        data = {'vertices': np.asarray([p.xy for p in vertices.values()]),
                'segments': np.asarray(segments)}

        t = triangle.triangulate(data, opts='p')  # triangulates

        # some polygons can not be triangulated, like those made up from a single segment
        # for these the variable t has no key 'triangles'
        if 'triangles' not in t.keys():
            return False

        # if the number of vertices has been changed, the dict vertices must be updated
        for v in t['vertices']:
            if Point(tuple(v)) not in vertices.values():
                vertices[len(vertices.keys())] = Point(tuple(v))

        for tri in t['triangles']:
            _poly = polygon_from_nodes([vertices[x].xy for x in tri])  # a polygon line from the nodes
            _poly.close()  # closing the line
            _tri = Triangle(_poly.segments)  # a triangle is formed from the polygons segments
            # checking validity
            if triangle_ok(polygon=self, triangle=_tri):
                # print('added')
                self._triangles.add(_tri)
            else:
                pass
                # print('not added')
                # print([_tri.is_internal_point(p=x) for x in self.point_set if x not in _tri.point_set])
                # self.plot(show=False)
                # _poly.plot(show=True)


        # # # this shows the result of the triangluation. non-valid triangles are possible
        # import triangle.plot as plot
        # ax = plt.axes()
        # plot.plot(ax, **t)
        # plt.show()
        #
        # # # this shows the resulting triangulation; non-valid triangles are thrown away
        # for _tri in self.triangles:
        #     _tri.plot(show=False)
        # _tri.plot(show=True)

    def plot_triangulation(self, algorithm=TRIANGULATION_ALGORITHM, force=False, show=False):
        # this plots the triangulation

        if not self.triangles or force:
            if algorithm is not None:
                if algorithm == 'delaunay':
                    self.triangulate_delaunay()
                elif algorithm == 'ear clipping':
                    self.ear_clipping()
                elif algorithm == 'triangle':
                    self.triangle_mesh()
                else:
                    raise Exception('Triangulation algorithm %s unknown' % algorithm)

            else:
                raise Exception('Please provide a name for the algorithm to be used')

        # as not all polygons can be triangulated
        if self.triangles:
            self.plot(show=False)
            for tri in self.triangles:
                tri.plot(show=False)
                tri.midpoint.plot(show=False, style='k.')
            tri.plot(show=show)

    def plot(self, segments=True, points=True, also_plot=(), show=False, annotate=True, title=None):
        if _plotting_available:
            if segments:
                [s.plot() for s in self.segments]
            if points:
                [p.plot() for p in self.point_set]
                self.segments[0].i.plot()
            if also_plot:
                for o in also_plot:
                    o.plot(show=False, annotate=annotate)

            # setting the axis
            plt.axis('tight')
            plt.axis('equal')
            if title is not None:
                plt.title(title)

            if show:
                plt.show()

        else:
            print('No plotting available')


class Triangle(PolygonLine):
    def __init__(self, args):
        super(Triangle, self).__init__(args)
        self.purge()
        self.order()
        if not self.is_triangle:
            raise Exception('non valid triangle')

    @property
    def detailled_area(self):
        """ the area calculated using different algorithms """
        A = self.area
        return {'delaunay': A, 'ear clipping': A, 'triangle': A}

    @property
    def area(self):
        """ area by the Heron-formula """
        s = sum([x.norm for x in self.segments]) / 2.
        _segs = self.segments
        discriminant = s * (s - _segs[0].norm) * (s - _segs[1].norm) * (s - _segs[2].norm)
        try:
            return math.sqrt(discriminant)
        except ValueError:  # if the discriminant would be zero
            if abs(discriminant) < EPS:
                return 0
            else:
                raise

    @property
    def quality(self):
        """ gives an idea on the perfectness of a triangle. lower is better. """
        return sum([abs(x - 60) for x in self.angles])

    @property
    def angles(self):
        """ internal angles of the triangle """
        _triangle = copy.deepcopy(self)
        _ret = []
        for p in itertools.combinations(_triangle.segments, 2):
            if p[0].i == p[1].j:
                p[1].reverse()
            elif p[1].i == p[0].j:
                p[0].reverse()
            elif p[1].j == p[0].j:  # endpoint is the same
                p[0].reverse()
                p[1].reverse()
            else:  # startpoints are the same
                pass
            _ret.append(p[0].angle_between_lines(other=p[1]))
        if not sum(_ret) - 180 < EPS:
            print(_ret)
        return _ret

    def congruent(self, other):
        """ all points are coincident, works vica versa """
        assert isinstance(self, Triangle) and isinstance(other, Triangle)

        if self.point_set == other.point_set:
            # print('congruent')
            return True

    def enclosed(self, other):
        """ other is in self, DOESN'T work vica versa """
        assert isinstance(self, Triangle) and isinstance(other, Triangle)

        # all points of other are internal for self
        if all(self.is_internal_point(x) for x in other.point_set):
            # print('other is fully enclosed in self')
            return True

    def overlap(self, other):
        """ partial or total overlap, works vica versa """
        assert isinstance(self, Triangle) and isinstance(other, Triangle)

        # some points of self are inside the other or vica versa
        if any(self.is_internal_point(x) for x in other.point_set) or \
                any(other.is_internal_point(x) for x in self.point_set):
            return True

        # there are intersecting lines but the intersection_point is not on an endpoint
        _common_parts = set()
        pairs = set()
        for ss in self.segments:
            for os in other.segments:
                if (ss, os) not in pairs and (os, ss) not in pairs:
                    pairs.add((ss, os))

        for pair in pairs:
            cp = pair[0].common_part(pair[1])
            cp_type = pair[0].common_part_type(pair[1])
            if cp_type == 'Point' and cp is not None and cp not in self.point_set + other.point_set:
                _common_parts.add(cp)
        if _common_parts:
            return True

        return False

    def touch(self, other):
        """ whatever is found in common_part, the points are all on the boundary """
        assert isinstance(self, Triangle) and isinstance(other, Triangle)

        pairs = set()
        for ss in self.segments:
            for os in other.segments:
                if (ss, os) not in pairs and (os, ss) not in pairs:
                    pairs.add((ss, os))

        _common_parts = set()
        for pair in pairs:
            cp = pair[0].common_part(pair[1])
            cp_type = pair[0].common_part_type(pair[1])
            if cp_type == 'LineSegment':
                _common_parts.add(cp)

        if len(_common_parts) == 1:
            return True
        else:
            return False

    def disjunct(self, other):
        """
        tells the polygonline self is disjunct from other. both self and other must be triangles.
        distjunct is, if:
        - not congruent
        - no partial overlap
        NOT disjunct in all other cases:
        - fully enclosed one way or another
        - they touch, that is have some sections that are at leas partially overlapping
        """
        assert isinstance(self, Triangle) and isinstance(other, Triangle)

        # total overlap
        if self.congruent(other):
            return False

        # partial overlap
        if self.overlap(other) or other.overlap(self):
            return False

        # fully enclosure
        if self.enclosed(other) and not other.enclosed(self):
            return False

        if self.touch(other):
            return False

        return True

    def tell_relpos(self, other):
        for posname in ['touch', 'overlap', 'disjunct', 'enclosed', 'congruent']:
            relpos = getattr(self, posname)(other)
            if relpos is True:
                print('%s: True' % posname)


class Delaunay2d:
    """
    this class creates a Delaunay triangulation of a PolygonLine.
    based on: http://code.activestate.com/recipes/579021-delaunay-triangulation/
    changes:
    - integration of the SR module
    - numpy as dependecy no longer required
    """

    def __init__(self, polygon):

        # my addition: to enable easier integration with the SR module
        self.polygon = polygon
        self.constraints = self.polygon.segments
        # self.nr_added_points = len(self.polygon.point_set)
        self.nr_added_points = 0

        # data structures
        self.points = copy.deepcopy(self.polygon.point_set)  # copy
        self.triangles = []  # cells
        self.edge2Triangles = {}  # edge to triangle(s) map
        self.boundaryEdges = set()
        self.boundaryEdges = set()
        self.appliedBoundaryEdges = None
        self.add_internal_points()
        self.points = list(set(self.points))

    def refine_boundaries(self):
        """ refines the boundaries to make triangulation more robust """
        _minlen = min([x.norm for x in self.constraints])
        for c in self.constraints:
            _divs = int(c.norm / _minlen) + 1
            for s in divide_line(c, _divs):
                self.points.append(s.i)
                self.points.append(s.j)

    def add_internal_points(self):
        """ adds internal points to make triangulation more robust """
        bb = self.polygon.bounding_box()
        _num_internal_point = 0
        while _num_internal_point < self.nr_added_points:
            ps = Point((random.uniform(min(bb[0]), max(bb[0])), random.uniform(min(bb[1]), max(bb[1]))))
            if self.polygon.is_internal_point(p=ps):
                self.points.append(ps)
                _num_internal_point += 1

    @property
    def all_segments(self):
        return itertools.chain.from_iterable([tri.segments for tri in self.polygon.triangles])

    @property
    def all_points(self):
        return itertools.chain.from_iterable([seg.ij for seg in self.all_segments])

    def update_points(self):
        """ makes sure that before triangulation the points are up-to-date """
        self.points = copy.deepcopy(self.polygon.point_set)  # copy
        self.add_internal_points()
        # self.refine_boundaries()
        self.points = list(set(self.points))
        # for p in self.points:
        #     p.plot(show=False)
        # p.plot(show=True)
        # exit()

    def triangulate(self):

        self.update_points()  # update points

        points = self.points
        # points = self.polygon.point_set
        # compute center of gravity
        cg = Point((0, 0))
        for pt in points:
            cg += pt
        cg._set_coords(coords=[x / len(points) for x in cg.xy])

        # sort
        def distanceSquare(pt):
            return distance(cg, pt) ** 2

        self.points.sort(key=distanceSquare)  # sorting the points by their distance from the centroid

        # create first triangle, make sure we're getting a non-zero detailled_area otherwise
        # drop the points
        area = 0.0
        index = 0
        stop = False
        while not stop and index + 2 < len(points):
            area = self.getArea(index, index + 1, index + 2)
            if abs(area) < EPS:
                del self.points[index]
            else:
                stop = True

        if index <= len(self.points) - 3:
            tri = [index, index + 1, index + 2]
            self.makeCounterClockwise(tri)
            self.triangles.append(tri)

            # boundary edges
            e01 = (tri[0], tri[1])  # edge between two points of the newly created triangle
            self.boundaryEdges.add(e01)
            e12 = (tri[1], tri[2])
            self.boundaryEdges.add(e12)
            e20 = (tri[2], tri[0])
            self.boundaryEdges.add(e20)

            e01 = self.makeKey(e01[0], e01[1])  # a tuple with the numbers in ascending order
            self.edge2Triangles[e01] = [0, ]
            e12 = self.makeKey(e12[0], e12[1])
            self.edge2Triangles[e12] = [0, ]

            e20 = self.makeKey(e20[0], e20[1])
            self.edge2Triangles[e20] = [0, ]

        else:
            # all the points fall on a line
            return

        # add additional points
        _pts = set(range(3, len(self.points)))
        for i in range(3, len(self.points)):
            self.addPoint(_pts.pop())

    def getTriangles(self):
        """
        @return triangles
        """
        return self.triangles

    def getEdges(self):
        """
        @return egdes
        """
        return self.edge2Triangles.keys()

    def getArea(self, ip0, ip1, ip2):
        """
        Compute the parallelipiped detailled_area
        @param ip0 index of first vertex
        @param ip1 index of second vertex
        @param ip2 index of third vertex
        """
        d1 = self.points[ip1] - self.points[ip0]
        d2 = self.points[ip2] - self.points[ip0]
        return d1.x * d2.y - d1.y * d2.x

    def isEdgeVisible(self, ip, edge):
        """
        Return true iff the point lies to its right when the edge points down
        @param ip point index
        @param edge (2 point indices with orientation)
        @return True if visible
        """
        area = self.getArea(ip, edge[0], edge[1])
        if area < EPS:
            return True
        return False

    def intersectsConstraint(self, edge):
        """
        return True if edge intersects a constraint defined for the triangulation
        """
        _edge = LineSegment(*[self.points[x] for x in edge])
        if any([x.intersection_on_segment(other=_edge) for x in self.constraints]):
            return True
        else:
            return False

    def makeCounterClockwise(self, ips):
        """
        Re-order nodes to ensure positive detailled_area (in-place operation)
        """
        area = self.getArea(ips[0], ips[1], ips[2])
        if area < -EPS:
            ip1, ip2 = ips[1], ips[2]
            # swap
            ips[1], ips[2] = ip2, ip1

    def flipOneEdge(self, edge):
        """
        Flip one edge then update the data structures
        @return set of edges to interate over at next iteration
        """

        # start with empty set
        res = set()

        # assume edge is sorted
        tris = self.edge2Triangles.get(edge, [])
        if len(tris) < 2:
            # nothing to do, just return
            return res

        try:
            iTri1, iTri2 = tris
        except ValueError:
            pp.pprint(tris)

            import time
            time.sleep(0.5)
            raise

        tri1 = self.triangles[iTri1]
        tri2 = self.triangles[iTri2]

        # find the opposite vertices, not part of the edge
        iOpposite1 = -1
        iOpposite2 = -1
        for i in range(3):
            if not tri1[i] in edge:
                iOpposite1 = tri1[i]
            if not tri2[i] in edge:
                iOpposite2 = tri2[i]

        # compute the 2 angles at the opposite vertices
        da1 = self.points[edge[0]] - self.points[iOpposite1]
        db1 = self.points[edge[1]] - self.points[iOpposite1]
        da2 = self.points[edge[0]] - self.points[iOpposite2]
        db2 = self.points[edge[1]] - self.points[iOpposite2]
        crossProd1 = self.getArea(iOpposite1, edge[0], edge[1])
        crossProd2 = self.getArea(iOpposite2, edge[1], edge[0])
        dotProd1 = self.dotproduct(da1, db1)
        dotProd2 = self.dotproduct(da2, db2)
        angle1 = abs(math.atan2(crossProd1, dotProd1))
        angle2 = abs(math.atan2(crossProd2, dotProd2))

        # Delaunay's test
        if angle1 + angle2 > math.pi * (1.0 + EPS):

            # flip the triangles
            #             / ^ \                    / b \
            # iOpposite1 + a|b + iOpposite2  =>   + - > +
            #             \   /                    \ a /

            newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
            newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

            # update the triangle data structure
            self.triangles[iTri1] = newTri1
            self.triangles[iTri2] = newTri2

            # now handle the topolgy of the edges
            # remove this edge
            del self.edge2Triangles[edge]

            # add new edge
            e = self.makeKey(iOpposite1, iOpposite2)
            self.edge2Triangles[e] = [iTri1, iTri2]

            # modify two edge entries which now connect to
            # a different triangle
            e = self.makeKey(iOpposite1, edge[1])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri1:
                    v[i] = iTri2
            res.add(e)

            e = self.makeKey(iOpposite2, edge[0])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri2:
                    v[i] = iTri1
            res.add(e)

            # these two edges might need to be flipped at the
            # next iteration
            res.add(self.makeKey(iOpposite1, edge[0]))
            res.add(self.makeKey(iOpposite2, edge[1]))

        return res

    def flipEdges(self):
        """
        Flip edges to statisfy Delaunay's criterion
        """

        # start with all the edges
        edgeSet = set(self.edge2Triangles.keys())

        continueFlipping = True

        while continueFlipping:

            #
            # iterate until there are no more edges to flip
            #

            # collect the edges to flip
            newEdgeSet = set()
            for edge in edgeSet:
                # union
                newEdgeSet |= self.flipOneEdge(edge)

            edgeSet = copy.copy(newEdgeSet)
            continueFlipping = (len(edgeSet) > 0)

    def addPoint(self, ip):
        """
        Add point to the triangulation
        @param ip point index
        """

        # collection for later updates
        boundaryEdgesToRemove = set()
        boundaryEdgesToAdd = set()

        for edge in self.boundaryEdges:

            if self.isEdgeVisible(ip, edge):
                # create new triangle
                newTri = [edge[0], edge[1], ip]
                newTri.sort()
                self.makeCounterClockwise(newTri)
                self.triangles.append(newTri)

                # update the edge to triangle map
                e = list(edge[:])
                e.sort()
                iTri = len(self.triangles) - 1
                self.edge2Triangles[tuple(e)].append(iTri)

                # add the two boundary edges
                e1 = [ip, edge[0]]
                e1.sort()
                e1 = tuple(e1)
                e2 = [edge[1], ip]
                e2.sort()
                e2 = tuple(e2)
                v1 = self.edge2Triangles.get(e1, [])
                v1.append(iTri)
                v2 = self.edge2Triangles.get(e2, [])
                v2.append(iTri)
                self.edge2Triangles[e1] = v1
                self.edge2Triangles[e2] = v2

                # keep track of the boundary edges to update
                boundaryEdgesToRemove.add(edge)
                boundaryEdgesToAdd.add((edge[0], ip))
                boundaryEdgesToAdd.add((ip, edge[1]))

        # update the boundary edges
        for bedge in boundaryEdgesToRemove:
            self.boundaryEdges.remove(bedge)
        for bedge in boundaryEdgesToAdd:
            bEdgeSorted = list(bedge)
            bEdgeSorted.sort()
            bEdgeSorted = tuple(bEdgeSorted)
            if len(self.edge2Triangles[bEdgeSorted]) == 1:
                # only add boundary edge if it does not appear
                # twice in different order
                self.boundaryEdges.add(bedge)

        # recursively flip edges
        flipped = True
        while flipped:
            flipped = self.flipEdges()

    def makeKey(self, i1, i2):
        """
        Make a tuple key such at i1 < i2
        """
        return tuple(sorted([i1, i2]))

    def dotproduct(self, a, b):
        """
        Dot product.
        """
        return a.x * b.x + a.y * b.y

    def show(self, width=500, height=500, showVertices=True, showContour=[]):

        from geometry import plt

        xmin = min([p.x for p in self.points])
        ymin = min([p.y for p in self.points])
        xmax = max([p.x for p in self.points])
        ymax = max([p.y for p in self.points])
        padding = 5
        w = width - 2 * padding
        h = height - 2 * padding

        for e in self.edge2Triangles:
            i1, i2 = e
            xp1 = padding + int(w * (self.points[i1].x - xmin) / (xmax - xmin))
            yp1 = padding + int(h * (ymax - self.points[i1].y) / (ymax - ymin))
            xp2 = padding + int(w * (self.points[i2].x - xmin) / (xmax - xmin))
            yp2 = padding + int(h * (ymax - self.points[i2].y) / (ymax - ymin))
            _ls = LineSegment(Point((xp1, yp1)), Point((xp2, yp2)))
            _ls.plot(show=False)

        if showVertices:
            for i in range(len(self.points)):
                xp = padding + int(w * (self.points[i].x - xmin) / (xmax - xmin))
                yp = padding + int(h * (ymax - self.points[i].y) / (ymax - ymin))
                _vp = Point((xp, yp))
                _vp.plot(show=False)

        if len(showContour) > 0:
            for i in range(len(showContour) - 1):
                xp1 = padding + int(w * (showContour[i].x - xmin) / (xmax - xmin))
                yp1 = padding + int(h * (ymax - showContour[i].y) / (ymax - ymin))
                xp2 = padding + int(w * (showContour[i + 1].x - xmin) / (xmax - xmin))
                yp2 = padding + int(h * (ymax - showContour[i + 1].y) / (ymax - ymin))
                _ls = LineSegment(Point((xp1, yp1)), Point((xp2, yp2)))
                _ls.plot(show=False)

        plt.axis('tight')
        plt.axis('equal')
        plt.show()


def segments_from_pointset(points):
    """ from a set of points it creates a list of LineSegments with the shortest possible length. """
    _minlength = 10e10

    for p in set(points):
        ps = [p]
        _rest = copy.deepcopy(points)
        for i in range(len(_rest)-1):
            _rest.remove(ps[-1])
            _rest = [x for x in points if x not in ps]
            ps.append(sorted(_rest, key=lambda x: distance(x, p))[0])

        ls = []
        for pindex, p in enumerate(ps[:-1]):
            _add = LineSegment(p, ps[pindex+1])
            if _add.isvalid:
                ls.append(_add)
        _sumlen = sum([x.norm for x in ls])
        if _sumlen < _minlength:
            _minlength = _sumlen
            _ret = ls

    return _ret

#
# def find_collinear_lines(polygon):
#     """
#     given a polygon this finds the lines that are made up from is_collinear points and returns them
#     the return value is a tuple. [0] is the line without the internal points, [1] is a list of lines
#     that are made by is_collinear points
#     """
#     _ret = []
#     for s in remove_collinear_points(polygon).segments:  # s are the LineSegments without the "internal" points
#         _internals = [x for x in polygon.point_set if s.is_internal_point(p=x)]  # the internal points
#         ls = segments_from_pointset(points=_internals)  # making a list of segments from the points
#         _ret.append((s, _internals, ls))  # return value: the line without internal points, internal points, new mlist of segments
#
#     return _ret
#
#
def remove_collinear_points(polygon):
    """
    removes all is_collinear points so that the endpoints of the lines previously having the is_collinear points remain
    thus is the polygon simplified
    """
    _c = copy.deepcopy(polygon)
    for p in [x.i for x in polygon.segments]:
        nb = _c.neighbour(p=p)
        _l = LineSegment(*nb)  # the line connecting the neighbours
        if _l.is_point_on_line(p=p):  # points are is_collinear, p is between the endpoints of _l
            _segs = _c.find_segment_by_point(point=p, whichend='both')
            _c._segments = [x for x in _c._segments if x not in _segs]
            _c._segments.append(_l)
            _c.order()
    return _c


# def dissect_polygon(polygon):
#     """ this dissects a polygon so that the resulting polynoms are triangulable """
#
#     _c = copy.deepcopy(polygon)
#     _ret = copy.deepcopy(polygon)
#     _simple = remove_collinear_points(_c)  # a simplified polygon: no is_collinear points
#     _c._triangles = set()
#     _ret._segments = []
#
#     _simple.triangulate()
#
#     # all segments of the triangulation
#     for tri in _simple._triangles:
#         for s in tri._segments:
#             for orig_line, points_in_line, segments_in_line in find_collinear_lines(_c):
#                 if set(s.ij) == set(orig_line.ij):
#                     _ret._segments += [w for w in tri._segments if set(w.ij) != set(s.ij)]
#                     _ret._segments += segments_in_line
#                     _ret.plot(show=True)
#
#     # pp.pprint(_ret)
#     _ret.purge()
#     _ret.plot(show=True)
#     exit()
#     return _ret


def rotationsmatrix(phi):
    """ rotation matrix. phi is degrees, positive phi rotates CW """
    phi = math.radians(phi)
    return [[math.cos(phi), 1 * math.sin(phi)], [-1 * math.sin(phi), math.cos(phi)]]


def simplify(polygon):
    """
    Given a polygon it removes outliers and their internal, external points
    """
    polygon._segments = remove_outlier_segments(polygon._segments)
    polygon._internal_point = [x for x in polygon._internal_point if x.segment in polygon.segments]
    polygon._external_point = [x for x in polygon._external_point if x.segment in polygon.segments]


def remove_outlier_segments(segs):
    """ removes from segs the segments that have Endpoints with one degree """
    _ch = PolygonLine(segs)
    _sp = [x for x in _ch.all_points if _ch.all_points.count(x) == 1]
    return [x for x in segs if not any([y in _sp for y in x.ij])]


def polygon_from_nodes(nodes, cw=False):
    """
    creates a polygon from a list of node coordinates. This is not a closed polygon, just a line
    """
    # list of nodes

    if cw:
        nodes = [nodes[0]] + nodes[::-1][:-1]  # reversed, keeping the first point

    ps = []
    for n in nodes:
        _add = Point(n)
        if _add.isvalid:
            ps.append(_add)

    # list of segemnts
    ls = []
    for pindex, p in enumerate(ps[:-1]):

        _add = LineSegment(p, ps[pindex+1])
        if _add.isvalid:
            ls.append(_add)

    return PolygonLine(ls)


def square(a=None, x0=None, y0=None, cw=False):
    """
    Square, side length a, position middle of the bottom line in x0, y0.
    cw or ccw on demand, result is invariant
    """
    ps = [(0, 0), (a/2., 0), (a/2., a), (-a/2., a), (-a/2., 0)]
    pl = polygon_from_nodes(nodes=ps, cw=cw)
    pl.close()
    pl.move(x0=x0, y0=y0)
    return pl


def arc(P=None, r=None, alpha_start=0, alpha_end=360, n=10, cw=False):
    """
    Circular arc
    :param P: Centerpoint, either list, tuple or Point
    :param r: radius
    :param alpha_start: angle to start
    :param alpha_end: angle to end
    :param n: number of SegmentLines
    :return: valid, continuous, purged PolygonLine object
    """
    assert 0 <= alpha_start < 360
    assert 0 < alpha_end <= 360
    assert alpha_start < alpha_end
    assert r > 0
    assert n > 0

    degstep = (alpha_end - alpha_start) / float(n)
    nodes = []

    # nodes = [[r * math.cos(math.radians((st * degstep) % 360)), r * math.sin(math.radians((st * degstep) % 360))] for st in range(n+1)]
    for st in range(n+1):
        _deg = (st * degstep) % 360
        nodes.append([r * math.cos(math.radians(_deg)), r * math.sin(math.radians(_deg))])
    pl = polygon_from_nodes(nodes=nodes, cw=cw)
    pl.rotate_about_origin(phi=alpha_start)
    if isinstance(P, Point):
        pl.move(x0=P.x, y0=P.y)
    elif isinstance(P, (type([]), type(()))):
        pl.move(x0=P[0], y0=P[1])
    else:
        raise Exception('Center of circular arc invalid, either list, tuple or Point object is excepted')

    return pl


def intersecting_lines(lines):
    """ provides a list where pairs of intersecting lines are listed """
    combs = itertools.combinations(lines, 2)
    _ret = [comb for comb in combs if comb[0].intersection_on_segment(comb[1])]
    return _ret


def convex_hull(point_set):
    """Computes the convex hull of a set of 2D point_set.

    Input: an iterable sequence of (x, y) pairs representing the point_set.
    Output: a PolyLine object with the points that form the convex hull of the input set,
    starting from the vertex with the lexicographically smallest coordinates.
    The resulting PolygonLine is convex, continuous, valid, ccw.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the point_set lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    point_set = sorted(set(point_set))

    # Boring case: no point_set or a single point, possibly repeated multiple times.
    if len(point_set) <= 1:
        return point_set

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the point_set are is_collinear.

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in point_set:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(point_set):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.

    ch = polygon_from_nodes(nodes=lower[:-1] + upper[:-1])
    ch.close()  # closes
    ch.purge()  # making sure points are identical, not equal
    return ch


def add(self, *other):
    """
    Adds open polygonlines to form a continuous, purged valid line
    polylines in other is added to self, that is the startpoint of other is attached to the endpoint
    of the previous (self, other etc.).
    If these points are equal, we get a longer polyline back.
    If they are not, they will be connected by a segment.
    If the PolygonLines have intersecting lines, an error is raised.

    """
    assert self.isvalid and all([x.isvalid for x in other])
    assert not (self.isclosed and all([x.isclosed for x in other]))

    _sum = self.segments
    for s in other:
        _sum += s.segments

    _ret = PolygonLine(_sum)
    _ret.purge()

    assert _ret.isvalid

    return _ret


def walk_segments(segments=(), startingpoint=None):
    """
    Starting at startingpoint walk a set of segments to see if they can make up a closed, ordered polygon.
    The polygon doesn't have to be ordered, but needs to be closed to walk it.
    """

    if startingpoint is None:
        startingpoint = segments[0].i

    _has_branch = False
    _not_visited = set(copy.deepcopy(segments))
    _aktsegment = [x for x in _not_visited if startingpoint == x.i][0]
    _visited = {_aktsegment}
    _aktpt = startingpoint
    while len(_not_visited) > 0:
        _nextsegment = [x for x in _not_visited if _aktpt == x.i and x != _aktsegment]  # lets hope for an ordered polygon
        if not _nextsegment:  # no ordered solution
            _nextsegment = [x for x in _not_visited if _aktpt in x.ij and x != _aktsegment]  # non-ordered solution
            if not _nextsegment:
                return False
        if len(_nextsegment) > 1:
            _has_branch = True
        _nextsegment = _nextsegment[0]

        _visited.add(_aktsegment)
        _not_visited.discard(_aktsegment)
        _aktsegment = copy.deepcopy(_nextsegment)
        _aktpt = [x for x in _aktsegment.ij if x != _aktpt][0]

        # check if we have a branchless closed loop?
        segs = copy.deepcopy(_visited)
        segs.add(_nextsegment)
        _poly = PolygonLine(segs)
        if _poly.order() and _poly.isclosed and not _has_branch:
            return _poly

    return False




def run():
    pass

if __name__ == '__main__':
    run()
