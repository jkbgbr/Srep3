# -*- coding: utf-8 -*-

from geometry.SR_point import EPS
from geometry.SR_polygon import PolygonLine, convex_hull, remove_outlier_segments, walk_segments
from cached_property import cached_property
import itertools
from geometry.SR_polygon import TRIANGULATION_ALGORITHM
from collections import Counter


class PolygonGroup(object):
    """
    Group of PolygonLine objects.
    Used to retreive information about the relative position of Polygons.
    Can be used to unite (approximately as in boolean unite) PolygonLines.
    """
    def __init__(self, args):
        if not args:
            args = (PolygonLine(None), PolygonLine(None))
        self._polygons = list(args)
        self._cycles = None
        self._graph = None

    @property
    def graph(self):
        raise NotImplementedError

    def __repr__(self):
        return 'PolygonGroup(%r)' % self._polygons

    @property
    def polygons(self):
        return self._polygons

    def polygon_from_repr(self, rep):
        _ret = [x for x in self.polygons if x == rep]
        if len(_ret) == 1:
            _ret = _ret[0]
        else:
            raise Exception('WTF')

        return _ret

    def name_from_repr(self, rep):
        """ retreives the names of the segments from the segment prints """

        def fill_vect(v):

            if isinstance(self, (Graph, BehalterGraph)):
                polygons = self.polygons
            elif isinstance(self, PolygonGroup):
                polygons = self.graph.polygons
            else:
                raise Exception('WTF')

            _vect = []
            for pl in v:
                _akt = [x for x in polygons if pl == x]
                if len(_akt) == 1:
                    _vect.append(_akt[0].name)
            return _vect

        if isinstance(rep, dict):
            _ret = {}
            for k, v in rep.items():
                _ret[k] = fill_vect(v)

        elif isinstance(rep, (list, set)):
            _ret = fill_vect(rep)

        elif isinstance(rep, PolygonLine):
            _ret = rep.name

        else:
            raise Exception('WTF')

        return _ret

    @property
    def cycles(self):
        """
        This returns the cycles of PolygonGroup using the numbering of the Polygongroup so that the segments
        are ordered to be continuous and provide the convex hull or the concave hull.
        """
        _ret = {}
        for k, v in self.segments_in_cycles.items():  # dict with cycle: points in cycle (integer numbered)
            segs = list(itertools.chain.from_iterable([x.wall_segments for x in v]))  # all segments in the cycle

            _ch = PolygonLine(remove_outlier_segments(segs))
            _ch.order()

            # _ret[k] = [self.points_numbered[x] for x in itertools.chain([x.i for x in _ch.segments])]  # Points, ordered
            # _ret[k] = list(itertools.chain([x.i for x in _ch.segments]))  # Points, ordered
            _ret[k] = _ch

        return _ret

    @property
    def disjunct_cycles(self):
        """ This returns the disjunct cycles of PolygonGroup using the numbering of the Polygongroup """
        _ret = {}
        for k, v in self.graph.disjunct_cycles.items():
            _ret[k] = [x.hauptteil for x in v]

        return _ret

    @property
    def overlaps_in_cycles_named(self):
        """ This returns the overlapping segments of PolygonGroup using the names of the polygons """
        _ret = {}
        for k, v in self.graph.overlaps_in_cycles.items():
            _mi = [x.hauptteil for x in v]
            _ret[k] = set()
            for m in _mi:
                _ret[k].add(self.name_from_repr(rep=m))
        return _ret

    @property
    def overlaps_in_cycles(self):
        """ This returns the overlapping segments of PolygonGroup using the numbering of the Polygongroup """
        _ret = {}
        for k, v in self.graph.overlaps_in_cycles.items():
            _ret[k] = [x.hauptteil for x in v]

        return _ret

    @property
    def contained_cycles(self):
        """ Since the numbers of the cycles is the same in the PoilygonGroup and the Graph... """
        return self.graph.contained_cycles

    @property
    def all_overlapping_segments(self):
        """ This returns all overlapping segments of PolygonGroup using the numbering of the Polygongroup """
        # self.graph.all_overlapping_segments  # cycles, using the numbering of the graph
        # todo
        return self.graph.all_overlapping_segments

    def segments_from_points(self, points):
        """ from a given set of points it returns the segments """
        # segments that contain all the points tht are found in the cycle
        _s = [x for x in self.all_segments if all([x.i in points and x.j in points])]

        # polygons of which at least one segment has all of its Points in the input points
        return [x for x in self.polygons if any([y in x.segments for y in _s])]

        # polygons that are entirely contained in the segments - doesn't work for split mantels
        # return [x for x in self.polygons if set(x.segments).issubset(set(_s))]

    @cached_property
    def segments_in_cycles(self):
        """ dict with cycle: segment list, using the segments of the PolygonGroup """
        sic = self.graph.segments_in_cycles
        _ret = {k: [x.hauptteil for x in v] for k, v in sic.items()}  # as object for use
        return _ret

    @property
    def isvalid(self):
        """
        Tells if the graph is valid
        nd is dict with the node degrees.
        nd = {0: 0, 1: 0, 2: n, 3: m, 4+: 0}.
        nd[0], nd[1]: not connected nodes, loose edges
        nd[2] >= 3, at least a triangle please
        nd[3] % 2 == 0, no odd number of branches
        nd[4+] == 0, max. 3 edges in a node
        """
        #
        _al = self.as_adjacency_list(as_numbered=True)
        nd = {}
        for n in range(10):
            nd[n] = sum([1 for x in _al.keys() if len(_al[x]) == n])

        if nd[0] > 0:  # has non-connected node
            print('has non-connected node')
            _ret = False
        elif nd[1] > 0:  # has loose edge node
            print('has loose edge')
            _ret = False
        elif nd[2] < 3:  # its not even a triangle
            print('its not even a triangle')
            _ret = False
        elif nd[3] % 2 == 1:  # has odd number of branches
            print('has odd number of branches')
            _ret = False
        elif sum([1 for x in range(4, 10) if nd[x] > 0]):  # has node with more than 3 edges connected
            print('has node with more than 3 edges connected')
            _ret = False
        else:
            _ret = True  # but there will be other chacks

        if not _ret:
            return _ret

        # all polygons, that are truly polygonlines must be valid
        _ret = all([x.isvalid for x in self.polygons])

        return _ret

    def plot(self, cycles=True, annotate=False):
        """ This plots all polygons in the group """

        if not cycles:
            self.polygons[0].plot(also_plot=self.polygons[1:], show=True, annotate=annotate)

        else:
            _cycles = self.cycles
            colors = ['red', 'green', 'blue', 'yellow']
            # self.polygons[0].plot(also_plot=self.polygons[1:], show=False, annotate=annotate)

            from geometry import plt, PG, PC

            ax = plt.gca()
            for k, v in _cycles.items():
                patches = []
                # first, the cycles - that is the wandung
                v.plot(show=False)
                # triangulating the raum
                v.triangulate()

                # plotting the triangles
                for tri in v.triangles:
                    tri.plot(show=False, also_plot=[tri.midpoint])
                    patches.append(PG([x.xy for x in tri.point_set], True))
                # creating the patches
                p = PC(patches, alpha=0.4, facecolors=colors[k], edgecolors=colors[k])
                ax.add_collection(p)
            plt.axis('tight')
            plt.axis('equal')
            plt.show()

    @property
    def all_segments(self):
        """ flattened list of segments """
        return itertools.chain.from_iterable(x.segments for x in self.polygons)
        # return sum([x.wall_segments for x in self.polygons], [])

    @property
    def representative_graph(self):
        """
        Creates a PolygonLine from a set of Polygonline segments.
        This PolygonLine is not necessarily valid, but can be used as a Graph of the PolygonGroup
        """
        # _sum = sum([x.segments for x in self.polygons], [])
        _sum = itertools.chain.from_iterable([x.segments for x in self.polygons])
        pl = PolygonLine(_sum)
        pl.purge()
        return pl

    # @property
    # def points_numbered(self):
    #     """
    #     assigns an integer value for all endpoints. No optimization in the numbering done.
    #     returns a dict, keys: points, values: integers
    #     """
    #     pl = self.representative_graph
    #
    #     return {k: v for k, v in zip(pl.point_set, range(len(pl.point_set)))}

    # @property
    # def numbers_pointed(self):
    #     """ the points_numbered in the other direction: given an int, it returns the Point """
    #     pl = self.representative_graph
    #     return {v: k for k, v in zip(pl.point_set, range(len(pl.point_set)))}

    @property
    def as_edgelist(self):
        """ Prepares a list of segments, where each node has a separate value (int), to be used in MCB """
        _conn = self.points_numbered
        # graph = sorted([(_conn[s.i], _conn[s.j]) for s in self.all_segments])
        graph = sorted([(s.i, s.j) for s in self.all_segments])
        return graph

    def as_adjacency_list(self, as_numbered=False):
        """ Prepares a dict with k: v = node: list(connected nodes) """
        pl = self.representative_graph

        if not as_numbered:
            _ret = {}
            for p in pl.point_set:
                _ret[p] = list(pl.neighbour(p))

        else:
            _ret = {}
            pn = self.points_numbered  # the dict for the mapping. rcm = None explicitely, and MUST be so.
            # as RCM needs a basic numbering
            for p in pl.point_set:
                _ret[pn[p]] = sorted([pn[x] for x in pl.neighbour(p)])

        # # # pre-sorting the values for each key to have the node degree in ascendeing order
        # for k, v in _ret.items():
        #     _ret[k] = sorted(_ret[k], key=lambda z: len(_ret[z]))

        return _ret

    def find_loops(self, graph=None, cycles=None, comblen=None):
        return None

    def branching_number(self, nodes):
        """ Tells how many in nodes have a degree more than two """
        return len([x for x in nodes if self.node_degree(_al=self.as_adjacency_list(), node=x)])

    def node_degree(self, _al, node):
        """ Tells the degree of node based on the adjacency list """
        assert node in _al.keys()
        return len(_al[node])

    # @property
    # def rcm(self):
    #     """
    #
    #     # rcm http://ciprian-zavoianu.blogspot.de/2009/01/project-bandwidth-reduction.html
    #     # the exampole is from this page
    #     _al = {
    #         1: (5,),
    #         2: (3, 6, 8),
    #         3: (2, 5),
    #         4: (7,),
    #         5: (1, 3),
    #         6: (2, 8),
    #         7: (4,),
    #         8: (2, 6),
    #     }
    #
    #     print([1,5,3,2,6,8,4,7]) # original
    #     print([0,4,2,1,5,7,3,6]) # zero-based result
    #     :return:
    #     """
    #
    #     # _al = {
    #     #     0: (4,),
    #     #     1: (2, 5, 7),
    #     #     2: (1, 4),
    #     #     3: (6,),
    #     #     4: (0, 2),
    #     #     5: (1, 7),
    #     #     6: (3,),
    #     #     7: (1, 5),
    #     # }
    #     _al = self.as_adjacency_list(as_numbered=True)
    #     # # pre-sorting the adjacency lists for each node to have the lists ordered by degree
    #     for k, v in _al.items():
    #         _al[k] = sorted(_al[k], key=lambda x: len(_al[x]))
    #
    #     r = []
    #     q = deque()
    #     _nodessorted = sorted(_al.keys(), key=lambda z: self.node_degree(_al, z))  # key numbers sorted based on degree
    #     for st1 in _nodessorted:
    #         if st1 not in r:
    #             r.append(st1)
    #         q.extend(_al[r[-1]])
    #         while q:
    #             c = q.popleft()
    #             if c not in r:
    #                 r.append(c)
    #             q.extend([x for x in _al[c] if x not in r])
    #
    #     return r


class BehalterGraph(PolygonGroup):
    def __init__(self, behalter=None):
        self.behalter = behalter
        super(BehalterGraph, self).__init__([p.graph for p in self.behalter.polygons])

    @property
    def cycles(self):

        if not self._cycles:  # calculate once

            _cycles = self.keep_base_polygons(polygons=self.find_loops())

            self._cycles = {k: v for k, v in zip(range(len(_cycles)), _cycles)}

        return self._cycles

    @property
    def disjunct_cycles(self):
        """
        This returns for each cycle the cycles that do not have acommon segment.
        If one cycle completely outside of another one, they are disjunct
        The returned value is a dict, where keys are the cycles, the values are lists of other cycles (if any)
        """
        _ret = {}

        if isinstance(self, (Graph, BehalterGraph)):
            _cycles = self.segments_in_cycles
        elif isinstance(self, PolygonGroup):
            _cycles = self.graph.segments_in_cycles
        else:
            raise Exception('WTF')

        for e in _cycles.keys():  # keys
            _ret[e] = []
            for z in [x for x in _cycles.keys() if x != e]:  # keys, except for "e"
                _intr = set(_cycles[e]).intersection(_cycles[z])  # intersection_point - common elements
                if len(set()) > 0:  # if the set is not empty
                    _ret[e].append(_intr)

        return _ret

    @property
    def contained_cycles(self):
        # """ for each cycle it gives back a list of cycles, that have all their points inside that cycle """
        _ret = {}
        if isinstance(self, (Graph, BehalterGraph)):
            _cycles = self.segments_in_cycles
        elif isinstance(self, PolygonGroup):
            _cycles = self.graph.segments_in_cycles
        else:
            raise Exception('WTF')

        for e in _cycles.keys():  # keys

            _ret[e] = []
            for z in [x for x in _cycles.keys() if x != e]:  # keys, except for "e"
                # _ch = itertools.chain.from_iterable([x.ij for x in _cycles[e]])
                _ch = sum([x.point_set for x in _cycles[e]], [])  # all points in cycle "e"
                _convhull = convex_hull([w.xy for w in _ch])  # as tuples for the convex hull
                # check if the points of z inside the convex hull of e are
                # _ch = itertools.chain.from_iterable([x.ij for x in _cycles[z]])
                _ch = sum([x.point_set for x in _cycles[z]], [])  # all points in cycle "z"
                # internal points for the convex hull
                _internals = [p for p in _ch if _convhull.is_internal_point(p)]

                if set(_internals) == set(_ch):  # totally contained
                    _ret[e].append(z)
                    # _convhull.plot(also_plot=_ch, show=True)  # plots the hull and the internal points
                elif not _internals:  # empty
                    pass
                # else:  # partial overlap
                #     pp.pprint(_internals)
                #     pp.pprint(_ch)
                #     raise Exception('WTF')

        return _ret
        #
        #
        #
        #
        # import time
        # time.sleep(0.5)




        # _ret = {}
        # if isinstance(self, (Graph, BehalterGraph)):
        #     _cycles = self.segments_in_cycles
        # elif isinstance(self, PolygonGroup):
        #     _cycles = self.graph.segments_in_cycles
        # else:
        #     raise Exception('WTF')
        #
        # for e in _cycles.keys():  # keys
        #
        #     _ret[e] = []
        #     for z in [x for x in _cycles.keys() if x != e]:  # keys, except for "e"
        #         # _ch = itertools.chain.from_iterable([x.ij for x in _cycles[e]])
        #         _ch = sum([x.point_set for x in _cycles[e]], [])  # all points in cycle "e"
        #         _convhull = convex_hull([w.xy for w in _ch])  # as tuples for the convex hull
        #         # check if the points of z inside the convex hull of e are
        #         # _ch = itertools.chain.from_iterable([x.ij for x in _cycles[z]])
        #         _ch = sum([x.point_set for x in _cycles[z]], [])  # all points in cycle "z"
        #         # internal points for the convex hull
        #         _internals = [p for p in _ch if _convhull.is_internal_point(p)]
        #
        #         if set(_internals) == set(_ch):  # totally contained
        #             _ret[e].append(z)
        #             # _convhull.plot(also_plot=_ch, show=True)  # plots the hull and the internal points
        #         elif not _internals:  # empty
        #             pass
        #         # else:  # partial overlap
        #         #     pp.pprint(_internals)
        #         #     pp.pprint(_ch)
        #         #     raise Exception('WTF')
        #
        # return _ret


    @property
    def overlaps_in_cycles(self):
        """
        This returns the overlaps in the cycles
        The returned value is a dict, where keys are tuples of the cycles and values the overlapping polygons
        """
        _ret = {}

        if isinstance(self, (Graph, BehalterGraph)):
            _cycles = self.segments_in_cycles
        elif isinstance(self, PolygonGroup):
            _cycles = self.graph.segments_in_cycles
        else:
            raise Exception('WTF')

        for e in _cycles.keys():  # keys
            for z in [x for x in _cycles.keys() if x != e]:  # keys, except for "e"
                _intr = list(set(_cycles[e]).intersection(_cycles[z]))  # intersection_point of sets - common elements
                # (z, e) would be the second time "e" and "z" are evaluated together
                if _intr and (z, e) not in _ret.keys():
                    _ret[tuple(sorted([e, z]))] = _intr  # sorted to make the order constant

        return _ret

    @property
    def all_overlapping_segments(self):
        """
        This returns a list that contains all segments that are overlapping parts.
        This list does not care which to cycles are the overlapping.
        """
        return sum(self.overlaps_in_cycles.values(), [])

    @cached_property
    def segments_in_cycles(self):
        """
        This gives back a list of segments, that are in the cycles found.
        Input is a list of cycles
        Then it finds the graphrepr_segments defined by these points
        """
        _ret = {}

        # for p in self.polygons:
        #     print(self.name_from_repr(p))
        #     print(p.segments)

        for k, v in self.cycles.items():
            #
            # print(k)
            # pp.pprint(v)
            # v.plot(show=True)
            # for p in self.polygons:
            #     p.plot(show=True, title=p.hauptteil.name)
            # p.plot(show=True)

            # finding the polygons that make up a segment
            # loose_subset: during finding the cycles some PolyGonlines are generated, these have the same segments as
            # the originals but may be slightly different (order, direction etc.)
            # any(): partial overlap, e.g. when one of the original polygons is split by a raumende
            _ret[k] = [x for x in self.polygons if v.loose_subset(other=x)
                       or any([y in x.segments for y in v.segments])]

        return _ret

    def find_loops(self, graph=None, cycles=None, comblen=None):
        """
        Finds the closed loops in self.polygons by trying all possible combinations with all possible lengths
        returned is a sorted list of closed, ordered, valid polygons that are potential candidates to be raum
        """

        # first loop: no input - initialization
        if graph is None:
            graph = list(set(self.all_segments))

        if cycles is None:
            cycles = []

        # first, try a walk-along starting at a random point. this may be faster if we have lots of segments but just
        # one cycle or two cycles.
        _polygon = walk_segments(segments=graph)
        while _polygon and _polygon not in cycles:  # a loop without branch
            cycles.append(_polygon)  # returns just a polygon
            graph = [x for x in graph if x not in _polygon.segments and x.line_reversed not in _polygon.segments]
            _polygon = walk_segments(segments=graph)

        # some loops may have been found
        comblen = 1
        while comblen <= len(list(graph)):
            comblen += 1
            # print('comblen: %d, nr. of combs: %d, graph size: %d' % (comblen, len(list(itertools.combinations(graph, r=comblen))), len(graph)))
            for p in itertools.combinations(graph, r=comblen):  # all combinations of comblen segments.

                # premature stop: all segment endpoints should be found two times
                # thus creating a Polygon can be
                # dict of values: nr. occurrencies of the segment endpoints
                c = Counter(itertools.chain.from_iterable(x.ij for x in p))
                if not set(c.values()) == {2}:  # should be 2 all the way
                    continue

                # make a polygonline from the current segments and order it
                _polygon = PolygonLine(list(p))
                # if it can be ordered then its is at least a closed loop
                if not _polygon.order():
                    print('cont: unorderable')
                    continue

                # accepted as cycles, if: valid, closed, previously not found
                if _polygon.isvalid and _polygon.isclosed and _polygon not in cycles:
                    cycles.append(_polygon)  # returns just a polygon

                    # check if some segments can be removed from graph.
                    # They can be removed, if the endpoints of the newly found polygon all are found only twice
                    # in the segments of graph, that is, only if the new polygon doesn't have a branching,
                    # can the segments be removed to reduce calculation time.
                    ps = itertools.chain.from_iterable(x.ij for x in _polygon.segments)  # endpoints of the polygon
                    branchpoints = (x for x in ps if list(itertools.chain.from_iterable(x.ij for x in graph)).count(x) != 2)
                    try:
                        next(branchpoints)
                    except StopIteration:  # a loop without branch was found, this can be removed.
                        graph = [x for x in graph if x not in _polygon.segments]
                        graph = [x for x in graph if x.line_reversed not in _polygon.segments]
                        continue

        polygons = sorted(cycles, key=lambda x: len(x.segments))

        return polygons

    def keep_base_polygons(self, polygons=None):
        """
        finds bases cycles by finding the set of triangulated polygons that are disjunct and have the maximal sum detailled_area
        """

        # if _removed is None:
        #     _removed = []

        # for p in polygons:
        #     p.plot(show=True)
        # p.plot(show=True)

        for p in polygons:
            p.triangulate(algorithm=TRIANGULATION_ALGORITHM)
            # p.plot_triangulation(show=True)
        # p.plot_triangulation(show=True)

        # shortcut for single-raum
        if len(polygons) == 1:
            return polygons

        elif len(polygons) > 1:

            _ret = set()

            # making sure no loosely equal polygons are kept
            _removed_polygons = set()
            for pe in polygons:
                for pz in [x for x in polygons if pe != x]:
                    if pe.loose_equal(pz) and pz not in _removed_polygons:
                        _removed_polygons.add(pe)
            polygons = [x for x in polygons if x not in _removed_polygons]

            # removing polygons that can be constructed from any number of other polygons
            # _removed_polygons = set()
            comblen = 1
            while comblen < len(polygons) + 1:
                comblen += 1
                combs = itertools.combinations(polygons, comblen)
                for comb in combs:
                    _all_tris = list(itertools.chain.from_iterable([x.triangles for x in comb]))
                    for p in polygons:
                        if abs(p.area - sum([tri.area for tri in _all_tris])) < EPS:
                            _removed_polygons.add(p)
                            polygons.remove(p)

            # creating pairs of polygons that overlap. later when making the combinations,
            # combinations that have any of these can be skipped right away
            _all_polygon_pairs = itertools.combinations(polygons, 2)
            _overlapping_polygonpairs = {x for x in _all_polygon_pairs if not x[0].disjunct(x[1])}

            for comblen in range(1, len(polygons)+1):

                # creating all combinations of polygons of given length, taking into account the overlapping pairs:
                # if any of the overlapping pairs is found in the combination, the combination can be omitted
                _all_combs = itertools.combinations(polygons, comblen)
                combs = (x for x in _all_combs if not any(set(y).issubset(set(x)) for y in _overlapping_polygonpairs))
                # print('comblen: %d, all: %d, reduced: %d, nr. of polygons: %d' % (comblen, len(_all_combs), len(combs), len(polygons)))

                for comb in combs:
                    # _ret.add(comb)  # the generator is empty - no overlapping triangles found

                    # all triangles from the combination of polygons
                    _all_tris = itertools.chain.from_iterable([x.triangles for x in comb])
                    pairs = (x for x in itertools.combinations(_all_tris, 2))

                    # if there are among the pairs some where the triangles somehow overlap, the polygons (cycles) in the
                    # comb overlap - and we don't need those
                    _congruents = (t for t in pairs if t[0].congruent(t[1]) or t[0].overlap(t[1]))
                    # checking if the generator is empty
                    try:
                        next(_congruents)  # there is something it the generator - pairs that are overlapping
                    except StopIteration:
                        _ret.add(comb)  # the generator is empty - no overlapping triangles found

            _ret = list(_ret)
            # sort the list by the number of polygons found, then by sum area of the triangles
            _ret.sort(key=lambda x: (len(x), sum([p.area for p in x])), reverse=True)

            # ordering the segments to provide a continuous polygon
            _rret = []
            for p in _ret[0]:
                p.order()
                _rret.append(p.closed_copy())

            # ordering the polygons, to be sure that equality correctly evaluated
            for p in polygons:
                p.order()

            # at this point the [0] element of _rret is the list of polygons that have the most possible polygons
            # and the larges possible area.
            # There are two cases.
            # #1 We have in [0] the inner polygons (e.g. multiple polygons inside)
            # #2 We have in [0] the outer polygon (e.g. simple doppelmantelbehalter)

            if len(polygons) != len(_rret):

                # the polygons that are not included in _rret; this is the outer shell or the inner shell
                _left_out = [x for x in polygons if x not in _rret and not any([x.loose_equal(y) for y in _rret])]
                if len(_left_out) != 1:
                    raise Exception('WTF')
                _left_out = _left_out[0]  # there is only one

                # _rret is a list of polygon-combinations.
                # for non-double-walled behalters _ret is a list, where the elements may not be overlapping
                # for the polygon combinations are previously filtered not to have such
                # the sorting provides a list with descending number of raums and sum area, so the first one
                # is the one we need.

                # for double walled behalters any of the elements may overlap with any others in general case
                # so here we should find out, which polygons are inside others
                # the criterium is:
                # - all of the endpoints are inside another polygon
                # we check all combinations; the order is important -> permutations, not combinations
                # from the combinations only those are interestig, where one of the polygons is already in _rret

                _to_add = []
                for oc in itertools.permutations(polygons, 2):  # all combinations of two
                    # oc[0].plot(show=False)
                    # oc[1].plot(show=True, title='%r' % [x.area for x in oc])
                    # print([x in _left_out for x in oc])
                    # if any([x in _left_out for x in oc]):  # that contain one from the polygons left out previously
                    if _left_out == oc[0]:  # that contain one from the polygons left out previously
                        if all((oc[0].is_internal_point(p=x) for x in oc[1].point_set)) and \
                                        oc[0] not in _to_add and oc[0] not in _rret:
                            # oc[0].plot(show=True, title='thisis, oc[0], %.2f' % oc[0].area)
                            _to_add.append(oc[0])
                        elif all((oc[0].is_external_point(p=x) for x in oc[1].point_set)) and \
                                        oc[0] not in _to_add and oc[0] not in _rret:
                            # oc[1].plot(show=True, title='thisis, oc[1], %.2f' % oc[1].area)
                            _to_add.append(oc[0])
                        else:
                            pass

                        # and any((x[0].overlap(x[1]) for x in itertools.combinations(_all_tris, 2))) \

                _rret += _to_add
                # _rret.insert(0, _to_add[0])
            # todo final check: no segments should be missing

            return _rret


class Graph(PolygonLine):
    def __init__(self, hauptteil=None):
        self.hauptteil = hauptteil
        self.name = 'Graph of %s' % self.hauptteil.name
        super(Graph, self).__init__(self.hauptteil.graphrepr_segments)

    def renew(self):
        """ re-creates the graph. To be used if the hauptteil has been changed (e.g. split) """
        self._segments = self.hauptteil.graphrepr_segments
        self._internal_point = []
        self._indifferent_point = []
        self._external_point = []
        self.add_intex_points()

    def add_intex_points(self):
        """
        Creates and adds internal or external points for the hauptteil of the Graph.
        These are generated based on the Polygon's existing Points, invards or outwards, using the hauptteil's
        intexpoint maker method.
        """

        for s in self.segments:
            for intex in self.hauptteil.intex_points[s].keys():
                for iep in self.hauptteil.intex_points[s][intex]:
                    if intex == 'internal':
                        self.add_internal_point(p=iep, segment=s)
                    elif intex == 'external':
                        self.add_external_point(p=iep, segment=s)
                    elif intex == 'indifferent':
                        self.add_indifferent_point(p=iep, segment=s)
                    else:
                        raise Exception('WTF')

        #     self.plot(show=False, also_plot=self.internal+self.external+self.indifferent+[self.midpoint])
        # self.plot(show=True, also_plot=self.internal+self.external+self.indifferent+[self.midpoint])
