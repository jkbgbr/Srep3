# -*- coding: utf-8 -*-

from geometry.SR_point import Point
from geometry.SR_line import LineSegment, point_on_line_at_ratio
from geometry.SR_polygon import PolygonLine, arc, intersecting_lines, simplify
from geometry.SR_polygon_group import PolygonGroup, BehalterGraph, Graph
import math
import pprint as pp
import itertools
import collections
import random
import copy
from cached_property import cached_property

KREMPE_SEGMENTS = 3
KALOTTE_SEGMENTS = 8
KEGELKREMPE_SEGMENTS = 3


# raum intex hivas eseten minden raum ketszer?
# todo: Union, s. dazu Handnotizen


def raum_dict_checker(dict_to_comapare=None, basedict=None):
    """
    this function checks if the values of a dict (intex_points) correspont to a given pattern
    the dict to be checked is 'dict_to_compare'
    the that is used as a reference is 'basedict'
    this has terrible code smell, but hey it works
    """

    def convert(dtc, _rs, p):
        for k, v in dtc.items():
            for kk, vv in dtc[k].items():
                _subst = set()
                for r in [x for x in _rs if x in vv]:
                    _subst.add(p[r])
                dtc[k][kk] = _subst
        return dtc

    # def dict_equal(d1, d2):
    #     return all([d1[k] == d2[k] for k in d1.keys()]) and all([d1[k] == d2[k] for k in d2.keys()])

    # unique values in the values()
    _raums = []
    for k, v in dict_to_comapare.items():
        for kk in dict_to_comapare[k].values():
            _raums.append(kk)
    _raums = list(set(itertools.chain.from_iterable(_raums)))

    # the 'independent' values used for comparison
    _iv = [chr(x) for x in range(97, min(123, 97+len(_raums)))]

    # the pairs to be checked. one of these must yield the same as the dict used for testing
    # a list of dicts is created for this
    pairs = [list(zip(x, _iv)) for x in itertools.permutations(_raums, len(_iv))]
    _repldict = []
    for pair in pairs:
        _repldict.append({})
        for p in pair:
            _repldict[-1][p[0]] = p[1]

    # converting the result dict once
    _basedict = copy.deepcopy(basedict)
    _basedict = convert(dtc=_basedict, _rs=_raums, p=_repldict[0])

    # converting and checking for all pairs
    _res = []
    _founds = []
    for pair in _repldict:
        _adict = copy.deepcopy(dict_to_comapare)
        _adict = convert(dtc=_adict, _rs=_raums, p=pair)

        # _res.append(dict_equal(d1=_adict, d2=_basedict))
        _res.append(_adict == _basedict)

        if _res[-1]:
            _founds.append(_adict)

    # the resulting list may have only one True value, all others must be False
    _trues = [x for x in _res if x]
    _falses = [x for x in _res if not x]
    # assertions
    if len(_trues) != 1:
        print('Problem: there should be exactly one solution to match the dicts, we have %d' % len(_trues))
        print(_res)
        print('The base dict')
        pp.pprint(_basedict)
        print('')
        print('The solutions found (if any)')
        for f in _founds:
            pp.pprint(f)
        if not _founds:
            pp.pprint(_adict)

        import time
        time.sleep(.2)
        raise Exception('WTF')

    if len(_trues) + len(_falses) != len(_res):
        print("Problem: sum of solutions doesn't match up: Trues: %d, Falses: %d, Sum: %d" %
              (len(_trues), len(_falses), len(_res)))
        print(_res)
        print('The base dict')
        pp.pprint(_basedict)
        print('')
        print('The solutions found (if any)')
        for f in _founds:
            pp.pprint(f)
        if not _founds:
            pp.pprint(_adict)

        import time
        time.sleep(.2)
        raise Exception('WTF')

    return True


def linenum():

    import fnmatch
    import os

    # dirnames = ['D:\\Srep2\\geometry\\',
    #             'D:\\Srep2\\SR_vessel\\',
    #             'D:\\Srep2\\tests\\']
    dirnames = ['C:\\Users\\X250-User\\PycharmProjects\\Srep2\\geometry\\',
                'C:\\Users\\X250-User\\PycharmProjects\\Srep2\\SR_vessel\\',
                'C:\\Users\\X250-User\\PycharmProjects\\Srep2\\tests\\']

    _ret = {}
    for dirn in dirnames:
        _dirsum = 0
        for root, dirname, filenames in os.walk(dirn):
            for filename in fnmatch.filter(filenames, '*.py'):
                with open(os.path.join(root, filename)) as f:
                    _c = 0
                    for i, l in enumerate(f):
                        if l not in ['\n', '\r\n'] and not l.startswith('#'):
                            _c += 1
                if filename not in _ret.keys():
                    _ret[filename] = _c
                else:
                    _ret[filename] += _c + 1
                _dirsum += _c
    return sum(_ret.values())


class Wall(LineSegment):
    def __init__(self, *args):
        super(Wall, self).__init__(*args)
        # self.typ = 'wall'
        self.plot_style = 'k-'


class NonWall(LineSegment):
    def __init__(self, *args):
        super(NonWall, self).__init__(*args)
        # self.typ = 'nonwall'
        self.plot_style = 'r:'


class Raum(object):
    def __init__(self, behalter=None, cycle=None):
        self.cycle = cycle  # nummer der cycle in behalter.graph.cycle
        self.behalter = behalter
        self._name = 'R%d' % (len(self.behalter.raume) + 1)
        self._aussere_wandung = []  # nested list: internal lists are cycles. Nested, b/c there may e.g. a shell side...
        self._innere_wandung = []  # nested list: internal lists are cycles. Multiple internal rooms are possibel

    def set_name(self, newname=None):
        self._name = newname

    @property
    def innere_wandung(self):
        """ The innere Wandung as hauptteile """
        return (x.hauptteil for x in itertools.chain.from_iterable(self._innere_wandung))

    @property
    def aussere_wandung(self):
        """ The aussere Wandung as hauptteile """
        return (x.hauptteil for x in itertools.chain.from_iterable(self._aussere_wandung))

    @property
    def ganze_wandung(self):
        """ The ganze Wandung as hauptteile """
        return itertools.chain.from_iterable([self.innere_wandung, self.aussere_wandung])

    @property
    def name(self):
        return self._name

    def graph_as_polygon(self):
        return PolygonLine(itertools.chain.from_iterable([x.graph.segments for x in itertools.chain(self.ganze_wandung)]))

    @property
    def as_polygon(self):
        """
        Creates a Polygon that is built based on the segments of all hauptteile
        in all the innere and aussere wandung. This is not necessarily a valid PolygonLine.
        """
        return PolygonLine(itertools.chain.from_iterable([x.wall_segments for x in itertools.chain(self.ganze_wandung)]))

    # def print_wandung_names(self, wandung):
    #     [print(self.behalter.name_from_repr(rep=x)) for x in wandung]


class Behalter(PolygonGroup):
    def __init__(self, args):
        super(Behalter, self).__init__(args)
        self.raume = []
        self._graph = BehalterGraph(behalter=self)
        self.cut_by_raumende()  # cuts by raumende, if any
        self.raum_erstellen()  # makes the rooms

    @property
    def graph(self):
        return self._graph

    def cut_by_raumende(self):
        """ splits the polygons by raumende """

        # splitting/cutting the unterlying polygons
        _re = [x for x in self.polygons if x.typ == 'raumende']
        if _re:
            for ht in self.polygons:
                for x in _re:
                    ht.split(ht=x, keep=False)
                # _res = [ht.split(ht=x, keep=False) for x in _re]  # a single T/F is returned

    def raum_erstellen(self):

        if self.raume:
            pass

        if not self.raume:

            # first we find the loops/base cycles/polygons that are
            # the raum, based on the Graph of the PolygonGroup
            # this is quite straightforward and works well by now

            # print('kezodnek a ciklusok')

            if isinstance(self, (Graph, BehalterGraph)):
                _cycles = self.behalter.cycles
            elif isinstance(self, PolygonGroup):
                _cycles = self.cycles
            else:
                raise Exception('WTF')

            # pp.pprint(_cycles)
            # print('megvannak a ciklusok')

            for c in _cycles:

                _r = Raum(cycle=[c], behalter=self)

                # find the segments/polygons that form the wandung of the raum
                # based on the segments included in the cycles (that is, segments of the Graph)
                # this is the connection between Graph and PolygonGroup
                if isinstance(self, (Graph, BehalterGraph)):
                    sic = self.segments_in_cycles
                elif isinstance(self, PolygonGroup):
                    sic = self.graph.segments_in_cycles
                else:
                    raise Exception('WTF')

                # aussere Wandung
                import itertools
                _r._aussere_wandung = list(itertools.chain([v for k, v in sic.items() if k == c]))  # the outer wall of the room is the cycle

                # for v in _r._aussere_wandung:
                #     # print(v)
                #     for vv in v:
                #         print(vv)
                #         vv.plot(show=False)
                #     vv.plot(show=True, title='aussere', annotate=True)

                # innere Wandung
                if self.contained_cycles[c]:

                    # print('there are internal walls')

                    _cycles_to_check = self.contained_cycles[c]
                    # first, check the cycles that are overlapping (if any)
                    _iw = []  # innere Wandung
                    _visited = []
                    for k, v in self.overlaps_in_cycles.items():  # k: tuple of cycles that have common segment, eg. (0, 1)
                        # first run
                        if not _visited:
                            _iw.append([])
                            for s in k:
                                _iw[-1] += [x for x in sic[s] if x not in self.all_overlapping_segments]

                        # later runs
                        for vindex, _v in enumerate(_visited):
                            if any([x in _v for x in k]):  # we saw this earlier!
                                for s in k:
                                    _iw[vindex] += [x for x in sic[s] if x not in self.all_overlapping_segments]
                            else:
                                _iw.append([])
                                for s in k:
                                    _iw[-1] += [x for x in sic[s] if x not in self.all_overlapping_segments]

                        _visited.append(k)
                        _cycles_to_check = [x for x in _cycles_to_check if x not in k]

                    _r.cycle.append(_iw)
                    _r._innere_wandung += _iw

                    # the cycles, that are not overlapping
                    for x in _cycles_to_check:
                        _r._innere_wandung.append(sic[x])

                else:
                    _r._innere_wandung = []

                # for v in _r._innere_wandung:
                #     # print(v)
                #     for vv in v:
                #         print(vv)
                #         vv.plot(show=False)
                #     vv.plot(show=True, title='innere')

                self.raume.append(_r)

        # for _r in self.raume:
        #     print(_r.name)
        #     print('innere')
        #     for s in _r.innere_wandung:
        #         s.plot(show=False)
        #         print(s)
        #     print('aussere')
        #     for s in _r.aussere_wandung:
        #         s.plot(show=False)
        #         print(s)
        #     s.plot(show=True)

    def graphical_testing(self):
        """ plots the vessel with random points that are blue if internal, red if external """
        ap = self.raume[0].as_polygon
        numpt = 1000
        randpts_x = (random.uniform(-10, 10) for x in range(numpt))
        randpts_y = (random.uniform(-10, 20) for x in range(numpt))
        ps = [Point((x, y)) for x, y in zip(randpts_x, randpts_y)]
        ap.plot(show=False)
        for p in ps:
            _int = ap.is_internal_point(plotit=False, p=p)
            if _int:
                p.plot(show=False, style='bo')
            else:
                p.plot(show=False, style='ro')
        ap.plot(show=True)

    def raum_intex(self):
        """
        Provides a dict with information if for a given hauptteil a given raum is internal or external
        Considered are only raums, that have the hauptteil in their internal or external boundary.
        A raum is considered internal, it a positive pressure in the raum works as internal pressure.
        A raum is considered external, if a positive pressure in the raum works as external pressure.
        Output example: [ma11: {'innenraum': R1, 'aussenraum': R2}, bo1: {'innenraum': R1, 'aussenraum': R2}] etc.
        """

        # making sure raume is constructed
        if not self.raume:
            self.raum_erstellen()

        _ret = {k: {'innenraum': set(), 'aussenraum': set(), 'indifferent': set()} for k in [x.name for x in self.polygons]}

        selfcopy = copy.deepcopy(self)

        # print('')
        for r in selfcopy.raume:
            # print('raum: %s' % r.name)
            raum_graph = r.graph_as_polygon()  # ein Polygon, representiert den Graph

            # here should made sure no outliers are present
            simplify(raum_graph)
            # raum_graph.plot(show=True, title=r.name)
            raum_graph.triangulate()
            # print('')
            for ht in (x for x in r.ganze_wandung if x in itertools.chain(r.aussere_wandung, r.innere_wandung)):
                # print(ht.name)
                ht_graph = ht.graph
                # print(ht_graph)
                # print('internal', ht_graph.internal)
                # print('external', ht_graph.external)
                # print('indifferent', ht_graph.indifferent)

                if ht_graph.internal:
                    # print('')
                    # print('vizsgalt: internal')
                    for s in ht_graph.segments:
                        # print(s)
                        _internal = [x for x in ht_graph.internal if x.segment == s]
                        if all((raum_graph.is_internal_point(p=P) for P in _internal)):
                            _ret[ht.name]['innenraum'].add(r.name)
                            # print(ht.name + '-hoz ' + r.name + ' mint internal')
                        # ht_graph.plot(show=True, also_plot=_internal)
                        # else:
                        #     print('NEM' + ht.name + '-hoz ' + r.name + ' mint internal')

                # pp.pprint(_ret)
                if ht_graph.external:
                    # print('')
                    # print('vizsgalt: external')
                    for s in ht_graph.segments:
                        # print(s)
                        _external = [x for x in ht_graph.external if x.segment == s]
                        if all((raum_graph.is_internal_point(p=P) for P in _external)):
                            _ret[ht.name]['aussenraum'].add(r.name)
                            # print(ht.name + '-hoz ' + r.name + ' mint external')
                        # else:
                        #     print('NEM' + ht.name + '-hoz ' + r.name + ' mint external')

                # pp.pprint(_ret)
                if ht_graph.indifferent:
                    # print('')
                    # print('vizsgalt: indifferent')
                    for s in ht_graph.segments:
                        _indifferent = [x for x in ht_graph.indifferent if x.segment == s]
                        if not all((raum_graph.is_internal_point(p=P) for P in _indifferent)):
                            _ret[ht.name]['indifferent'].add(r.name)
                            # print(ht.name + '-hoz ' + r.name + ' mint indifferent')
                        # else:
                        #     print('NEM' + ht.name + '-hoz ' + r.name + ' mint indifferent')

        return _ret

    def plot(self, cycles=True, annotate=False):
        """ This plots a Behalter, based on the raums """
        from geometry import plt, PG, PC
        for rindex, r in enumerate(self.raume):
            _segs = [w.wall_segments for w in itertools.chain(r.innere_wandung, r.aussere_wandung)]
            # for w in itertools.chain(r.innere_wandung, r.aussere_wandung):
            #     _segs.append(list(itertools.chain(w.wall_segments)))
            _poly = PolygonLine(itertools.chain.from_iterable(_segs))

            colors = ['red', 'green', 'blue', 'yellow', 'red', 'green', 'blue', 'yellow', 'red', 'green', 'blue', 'yellow', 'red', 'green', 'blue', 'yellow']

            ax = plt.gca()
            patches = []
            # first, the cycles - that is the wandung
            if cycles:
                _poly.plot(show=False)
            # triangulating the raum
            _poly.triangulate()

            # plotting the triangles
            for tri in _poly.triangles:
                # tri.plot(show=False, also_plot=[tri.midpoint])
                patches.append(PG([x.xy for x in tri.point_set], True))
            # creating the patches
            p = PC(patches, alpha=0.5, facecolors=colors[rindex])# , edgecolors=colors[rindex]
            ax.add_collection(p)
        plt.axis('tight')
        plt.axis('equal')
        plt.show()


class Hauptteil(PolygonLine):
    def __init__(self, typ=None, segments=()):
        super(Hauptteil, self).__init__(segments)
        self.typ = typ
        self._graph = None
        self.purge()

    @classmethod
    def from_dict(cls, adict):
        raise NotImplementedError

    @property
    def intex_points(self):
        """ internal, external, indifferent Point maker """
        raise NotImplementedError

    @property
    def characteristic_size(self):
        """
        The size characteristic for the Hauptteil
        Torispherical, Flachbode: the diameter
        Mantel, Kegelmantel: smallest diameter
        """
        raise NotImplementedError

    @property
    def graph(self):
        return self._graph

    @property
    def isvalid(self):
        """
        A Hauptteil is valid if it is valid an its graph is valid.
        Since some subclasses are not valid polylines (e.g Flachboden) the evaluation is individual.
        """
        raise NotImplementedError

    @property
    def graphrepr_segments(self):
        return self.wall_segments

    @property
    def wall_segments(self):
        return [x for x in self.segments if x.__class__.__name__ == 'Wall']
        # return [x for x in self.segments if x.typ == 'wall']

    @property
    def nonwall_segments(self):
        # return [x for x in self.segments if x.typ == 'nonwall']
        return [x for x in self.segments if x.__class__.__name__ == 'NonWall']

    @property
    def polygonline(self):
        """ the Hauptteil object as Polygonline """
        return PolygonLine(self.segments)

    def split(self, v=None, ht=None, keep=True):
        """
        splits the Polygonline at given vertical position or with given segment
        This is like a horizointal cut. Results a new graph but the cycles etc. stay the same.
        A horizontal line (e.g. Flachboden) is not cut when given a cut height same as its v position.
        If the cutting line touches another line, nothing happens
        v: cut at this vertical position
        seg: use this segmentline to cut. the semgentline must have on-segment intersection_point points
        keep: T/F, works only if a ht is given for the cut.
        if True, the cutting line is kept as is, False: parts that are internal (of the cut ht) will be
        thrown away and the ht gets redefined.
        """

        assert not ht == v is None

        # no cuts in TorispherischesBoden
        if 'TorispherischesBoden' in [x.__name__ for x in self.__class__.__bases__]:
            return False

        if ht is not None:
            assert isinstance(ht, RaumEnde)  # if we cut, cut only with RaumEnde
            _cl = ht.graph.segments[0]
            _keep = keep

        elif v is not None:
            _keep = True
            _bb = self.graph.bounding_box()  # bounding box
            _cl = LineSegment(Point((min(_bb[0]) - 0.5 * abs(max(_bb[0])-min(_bb[0])), v)),
                              Point((max(_bb[0]) + 0.5 * abs(max(_bb[0])-min(_bb[0])), v)))  # cutting line

        else:
            print(ht)
            print(v)
            import time
            time.sleep(1)
            raise Exception('WTF')

        cut_lines = intersecting_lines([_cl] + self.graph.segments)
        # cut lines is a nested list
        # each element is a tuple with: (cutter segment, cut segment)

        # case no intersecting lines
        if not cut_lines:
            return False

        # removed are the cuts that have intersection_point points at the endpoints (not internal points)
        cut_lines = [x for x in cut_lines if x[0].intersection_on_segment(x[1]) and not x[0].intersection_on_endpoint(x[1])]
        # cut_lines = [x for x in cut_lines if not any([x[0].is_internal_point(z) for z in x[1].ij])]
        if not cut_lines:
            return False

        # case there are some lines cut
        _new_segments = []
        gc = copy.deepcopy(self.graph)  # graph copy
        for sl in gc._segments:
            if not [x for x in list(itertools.chain.from_iterable(cut_lines)) if x == sl]:  # sl is not cut by _cl
                _new_segments.append(sl)
            else:
                # the segment can be found in cut_lines: the segment gets removed and redefined
                for _i in [x for x in cut_lines if sl in x]:
                    _ip = _i[0].intersection_point(_i[1])  # intersection_point Point

                    # remove the segment that has been split
                    if sl in self._segments:
                        self._segments.remove(sl)  # remove the segment that has been split

                    # new segments, matching the old type
                    if sl.__class__.__name__ == 'Wall':
                        self._segments.append(Wall(Point(sl.i.xy), Point(_ip.xy)))  # new segment one
                        self._segments.append(Wall(Point(_ip.xy), Point(sl.j.xy)))  # new segment two
                    elif sl.__class__.__name__ == 'NonWall':
                        self._segments.append(NonWall(Point(sl.i.xy), Point(_ip.xy)))  # new segment one
                        self._segments.append(NonWall(Point(_ip.xy), Point(sl.j.xy)))  # new segment two
                    else:
                        raise Exception('What else could it be')

        # making sure the numbering of the graphs is renewed
        # without this the node numbering is screwed up
        self.graph.renew()

        if not _keep:
            # the cutter segment will be redefined: parts, that are internal segments of the cut ht will be thrown away.
            _pts = list(_cl.ij)  # endpoints of the cutter line
            for _i in [x for x in cut_lines if _cl in x]:
                _pts.insert(1, _i[0].intersection_point(_i[1]))  # intersection_point Points are inserted

            # currently the only possibility is to have a ring-formed raumende.
            # the two outer parts are kept as Wall, the internal as NonWall
            ht._segments = [Wall(_pts[0], _pts[1]), NonWall(_pts[1], _pts[2]), Wall(_pts[2], _pts[3])]

            # making sure the numbering of the graphs is renewed
            # without this the node numbering is screwed up
            ht.graph.renew()

        return True


class Mantel(Hauptteil):
    def __init__(self, name=None, pos_v=None, pos_h=0, d=(), h=None):
        # super(Mantel, self).__init__(typ='mantel')
        self.name = name
        self.pos_v = pos_v
        self.pos_h = pos_h
        self.d = d
        self.h = h
        self._segments = self.make_segments
        super(Mantel, self).__init__(typ='mantel', segments=self.segments)
        self._graph = Graph(hauptteil=self)
        self.graph.add_intex_points()

    @classmethod
    def from_dict(cls, adict):
        try:
            mantel = Mantel(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
                h=adict['h'],
            )
            assert all([hasattr(mantel, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Mantel: %s' % e)

        return mantel

    @property
    def intex_points(self):
        """
        internal, external, indifferent Point maker for the Graph
        the Points are defined in a local coordinate system, which is:
        x - horizontal, y - vertical. origo: Axis.
        three internal points and one external.
        """

        _res = {}

        for s in self._graph.segments:
            _ret = []
            _res[s] = {'internal': [], 'external': [], 'indifferent': []}
            # the endpoints of a segment in ascending order vertically; lower, upper points
            _lp, _up = sorted(s.ij, key=lambda x: x.y)
            # "unit" for _lp, _up
            _e = [abs(min((_up.y - _lp.y) / 10., min((_p.x - self.pos_h) / 50000., .0020))) for _p in (_lp, _up)]
            # we create four points, two at each end of the segments
            for ratio in [0.01, 0.99]:
                _tp = point_on_line_at_ratio(line=s, ratio=ratio)
                _ret.append(Point((_tp.x + _e[0], _tp.y)))
                _ret.append(Point((_tp.x - _e[0], _tp.y)))
                # _tp.plot(show=False, style='k.')

            # _ret.append(Point((_lp.x + _e[0], _lp.y + 2 * _e[1])))
            # _ret.append(Point((_lp.x - _e[0], _lp.y + 2 * _e[1])))
            # _ret.append(Point((_up.x + _e[0], _up.y - 2 * _e[1])))
            # _ret.append(Point((_up.x - _e[0], _up.y - 2 * _e[1])))
            #
            # s.plot(show=False)
            # _lp.plot(show=False)
            # _up.plot(show=False)
            # for r in _ret:
            #     r.plot(show=False)
            # r.plot(show=True)

            for r in _ret:
                if self.is_internal_point(p=r) and not self.is_external_point(p=r):
                    _res[s]['internal'].append(r)
                elif not self.is_internal_point(p=r) and self.is_external_point(p=r):
                    _res[s]['external'].append(r)
                else:
                    raise Exception('Telling if Point internal or external impossible')

        #     self.plot(show=False, also_plot=_res[s]['internal'] + _res[s]['external'])
        # self.plot(show=True, also_plot=_res[s]['internal'] + _res[s]['external'])
        return _res

    @property
    def characteristic_size(self):
        return min(self.d) / 2.

    @property
    def isvalid(self):

        # checking if the graph is valid
        _s = self.graph.segments  # shorthand
        # two segments
        if len(_s) != 2:
            return False

        # coordinates
        _ok = []
        _ok.append((_s[0].i.y == _s[1].j.y))  # lower points y coordinates
        _ok.append((_s[0].j.y == _s[1].i.y))  # upper points y coordinates
        _ok.append((_s[0].i.x + _s[1].j.x) / 2 == (_s[0].j.x + _s[1].i.x) / 2)  # axis vertical
        if not all(_ok):
            return False

        # checking if the PolygonLine made up of the segments is valid
        return self.polygonline.isvalid

    @property
    def make_segments(self):
        p1 = Point((self.pos_h + self.d[0] / 2., self.pos_v))
        p2 = Point((self.pos_h + self.d[1] / 2., self.pos_v + self.h))
        p3 = Point((self.pos_h - self.d[1] / 2., self.pos_v + self.h))
        p4 = Point((self.pos_h - self.d[0] / 2., self.pos_v))
        ri_vert = Wall(p1, p2)
        le_vert = Wall(p3, p4)
        lo_hor = NonWall(p4, p1)
        hi_hor = NonWall(p2, p3)
        return [ri_vert, hi_hor, le_vert, lo_hor]


class TorispherischesBoden(Hauptteil):
    def __init__(self, name=None, pos_v=None, pos_h=0., d=None, h1=None, orient=1):
        self.name = name
        self.pos_v = pos_v
        self.pos_h = pos_h
        self.d = d
        self.h1 = h1
        self.orient = orient
        self._segments = self.make_segments
        super(TorispherischesBoden, self).__init__(typ='boden', segments=self.segments)
        self._graph = Graph(hauptteil=self)
        self.graph.add_intex_points()

    @property
    def intex_points(self):
        """
        internal, external, indifferent Point maker for the Graph
        the Points are defined in a local coordinate system, which is:
        x - horizontal, y - vertical. origo: Axis.
        three internal points and one external.
        """

        _res = {}
        for s in self._graph.segments:
            _res[s] = {'internal': [], 'external': [], 'indifferent': []}
            _pts = sorted(s.ij, key=lambda x: x.x)
            if min([x.x for x in _pts]) < self.pos_h:
                _vorient = -1
            elif max([x.x for x in _pts]) > self.pos_h:
                _vorient = 1
            else:
                raise Exception('WTF')
            _res[s]['internal'].append(Point((self.pos_h + _vorient * 0.45 * self.d, self.pos_v + self.orient * 0.05 * self.h_graph)))
            _res[s]['internal'].append(Point((self.pos_h + _vorient * 0.05 * self.d, self.pos_v + self.orient * 0.85 * self.h_graph)))
            _res[s]['external'].append(Point((self.pos_h + _vorient * 0.25 * self.d, self.pos_v + self.orient * 0.55 * self.h_graph)))

        #     self.graph.plot(show=False, also_plot=_res[s]['internal'] + _res[s]['external'])
        # self.graph.plot(show=True, also_plot=_res[s]['internal'] + _res[s]['external'])

        return _res

    @property
    def characteristic_size(self):
        return self.d / 2.

    @property
    def h_graph(self):
        return max(0.1, self.d / 1000.)

    @property
    def graphrepr_segments(self):
        return [Wall(Point((-self.d/2. + self.pos_h, self.pos_v)), Point((self.pos_h, self.pos_v + self.orient * self.h_graph))),
                Wall(Point((self.pos_h, self.pos_v + self.orient * self.h_graph)), Point((self.d / 2. + self.pos_h, self.pos_v)))]

    @property
    def isvalid(self):
        return self.polygonline.isvalid

    @property
    def make_segments(self):

        d = self.d
        pos_h = self.pos_h
        pos_v = self.pos_v
        h1 = self.h1

        if self.__class__.__name__ == 'Klopperboden':
            R1 = 1. * d  # kalotte
            R2 = 0.1 * d  # krempe
        elif self.__class__.__name__ == 'Korbbogenboden':
            R1 = 0.8 * d  # kalotte
            R2 = 0.154 * d  # krempe
        else:
            raise Exception('WTF')

        h = math.sqrt(((R1 - R2) ** 2) - ((d / 2. - R2) ** 2))  # Abstand TL bis Mittelpkt. Kalottenkugel
        alpha = math.degrees(math.acos(h/(R1 - R2)))  # halbe öffnungswinkel Kalottenkugel
        delta = 90 - alpha  # Öffnungswinkel Krempe

        # kalotte
        _mp_R1 = (0, h1 - h)
        # creates a PolygonLine
        _kal = arc(P=Point(_mp_R1), r=R1, alpha_start=90-alpha, alpha_end=90+alpha, n=KALOTTE_SEGMENTS)

        # krempe rechte Seite

        import copy
        _mp_R2 = (copy.deepcopy(_mp_R1)[0] + (R1 - R2) * math.sin(math.radians(alpha)),
                  copy.deepcopy(_mp_R1)[1] + (R1 - R2) * math.cos(math.radians(alpha)))
        _kre_re = arc(P=Point(_mp_R2), r=R2, alpha_start=0, alpha_end=delta, n=KREMPE_SEGMENTS)
        # krempe linke Seite
        _kre_li = arc(P=Point((-_mp_R2[0], _mp_R2[1])), r=R2, alpha_start=90 + alpha, alpha_end=180, n=KREMPE_SEGMENTS)

        zylboard_re = PolygonLine([Wall(Point((d/2., 0)), Point((d/2., h1)))])
        krempe_re = PolygonLine([Wall(x.i, x.j) for x in _kre_re.segments])
        kalotte = PolygonLine([Wall(x.i, x.j) for x in _kal.segments])
        krempe_li = PolygonLine([Wall(x.i, x.j) for x in _kre_li.segments])
        zylboard_li = PolygonLine([Wall(Point((-d/2., h1)), Point((-d/2., 0)))])
        RN = NonWall(Point((-d/2., 0)), Point((d/2., 0)))

        if self.orient == -1:
            zylboard_li.rotate_about_origin(phi=180)
            zylboard_re.rotate_about_origin(phi=180)
            krempe_li.rotate_about_origin(phi=180)
            krempe_re.rotate_about_origin(phi=180)
            kalotte.rotate_about_origin(phi=180)
            RN.reverse()

        zylboard_li.move(x0=pos_h, y0=pos_v)
        zylboard_re.move(x0=pos_h, y0=pos_v)
        krempe_li.move(x0=pos_h, y0=pos_v)
        krempe_re.move(x0=pos_h, y0=pos_v)
        kalotte.move(x0=pos_h, y0=pos_v)
        RN.move(x0=pos_h, y0=pos_v)

        _ret = zylboard_re.segments + krempe_re.segments + \
               kalotte.segments + krempe_li.segments + zylboard_li.segments + [RN]

        return _ret


class Klopperboden(TorispherischesBoden):
    def __init__(self, name=None, pos_v=None, pos_h=0, d=None, h1=None, orient=1):
        super(Klopperboden, self).__init__(name=name, pos_v=pos_v, pos_h=pos_h, d=d, h1=h1, orient=orient)

    @classmethod
    def from_dict(cls, adict):
        try:
            klopperboden = Klopperboden(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
                h1=adict['h1'],
                orient=adict['orient'],
            )
            assert all([hasattr(klopperboden, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Klopperboden: %s' % e)

        return klopperboden


class Korbbogenboden(TorispherischesBoden):
    def __init__(self, name=None, pos_v=None, pos_h=0., d=None, h1=None, orient=1):
        super(Korbbogenboden, self).__init__(name=name, pos_v=pos_v, pos_h=pos_h, d=d, h1=h1, orient=orient)

    @classmethod
    def from_dict(cls, adict):
        try:
            korbbogenboden = Korbbogenboden(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
                h1=adict['h1'],
                orient=adict['orient'],
            )
            assert all([hasattr(korbbogenboden, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Korbbogenboden: %s' % e)

        return korbbogenboden


class Flachboden(Hauptteil):
    def __init__(self, name=None, pos_v=None, d=None, pos_h=None):
        self.name = name
        self.d = d
        self.pos_h = pos_h
        self.pos_v = pos_v
        fb = Wall(Point((-self.d/2., 0)), Point((self.d/2., 0)))
        fb.move(x0=pos_h, y0=pos_v)
        super(Flachboden, self).__init__(typ='flachboden', segments=[fb])
        self._graph = Graph(hauptteil=self)
        self.graph.add_intex_points()

    @classmethod
    def from_dict(cls, adict):
        try:
            flachboden = Flachboden(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
            )
            assert all([hasattr(flachboden, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Flachboden: %s' % e)
        return flachboden

    @property
    def intex_points(self):
        """
        internal, external, indifferent Point maker for the Graph
        the Points are defined in a local coordinate system, which is:
        x - horizontal, y - vertical. origo: Axis.
        two indifferent points
        """
        _ret = {}
        for s in self.wall_segments:
            _ret[s] = {'internal': [], 'external': [], 'indifferent': []}
            # indifferent point, above
            _ret[s]['indifferent'].append(Point((self.pos_h, 0.01 * self.d + self.pos_v)))
            # indifferent point, below
            _ret[s]['indifferent'].append(Point((self.pos_h, -0.01 * self.d + self.pos_v)))
        return _ret

    @property
    def isvalid(self):
        """ valid if: valid as segment """
        _seg = self.segments[0]
        return _seg.isvalid

    @property
    def characteristic_size(self):
        return self.d / 2.

    def is_internal_point(self, p=Point()):
        """ tells if given point is internal. True if internal """
        assert self.isvalid
        assert p.isvalid

        # if the point is one of the Points of the polygon - False.
        if p in self.point_set:
            return False

        _seg = self.segments[0]
        if _seg.is_internal_point(p=p):  # called method of LineSegment object, self is a PolygonLine instance
            return True  # internal points are Points on the line. we should probably never get here
        else:
            return False


class Rohrboden(Flachboden):
    def __init__(self, name=None, pos_v=None, d=None, pos_h=None):
        super(Rohrboden).__init__(name=name, pos_v=pos_v, d=d, pos_h=pos_h)
        self.graph.add_intex_points()

    @classmethod
    def from_dict(cls, adict):
        try:
            rohrboden = Rohrboden(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
            )
            assert all([hasattr(rohrboden, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Rohrboden: %s' % e)

        return rohrboden

    def intex_points(self):
        raise NotImplementedError


class KegelKrempe(Hauptteil):
    def __init__(self, name=None, h_non_kegelside=None, h_kegelside=None, r=0., hts=None, kegel=None):
        """ 
        This creates a KegelKrempe at a junction of a Mantel and a Kegel. only converging are supported for now.
        Provide in hts the Mantel above or the Mantel above and under as a tuple.
        h_kegelside is the length of the abklinglange on the kegel
        h_non_kegelside is the length of the abklinglange on the ht connecting to the kegel
        hts[0] is a must, it must be an instance of Mantel
        hts[1] is optional, it must be either Mantel or StandZarge
        kegel is the kegel, it is a must
        Kegel is the cone connecting.
        
        Creates a Krempe.
        h_oben: Abstand der oberen RN von der Kegel-Mantel theoretische RN.
        h_unten: Abstand der unteren RN von der Kegel-Mantel theoretische RN.
        für r>0 sind beide erst dann wirksam, wenn die Längen grösser sind als die Abstände die sich aus Da und r ergeben.
        falls diese kleiner sind, werden die Endpunkte der Kurve verwendet
        für r=0 werden wirksam, falls die beide grösser als null sind. Falls diese gleich null sind, sollte
        keine Krempe definiert verwendet werden.
        
        """

        # todo: divergierende kegel
        # todo: kegel with cone "up"

        self.name = name
        self.h_oben = h_non_kegelside  # vertical distance of the RN in the upper ht
        self.h_unten = h_kegelside  # vertical distance of the RN in the lower ht
        self.r = r  # radius der Krempe
        self.hts = sorted(hts, key=lambda x: x.pos_v)
        self.kegel = kegel
        self.assert_checking()
        # todo: should say if ht_oben and ht_unten are not mantel
        self._graphrep = []  # will be defined in make_segments
        self._segments = self.make_segments
        super(KegelKrempe, self).__init__(typ='kegelkrempe', segments=self.segments)
        self._graph = Graph(hauptteil=self)
        self.graph.add_intex_points()
        self.order()
        self.purge()
        self.graph.renew()

    @property
    def ht_oben(self):
        return self.hts[-1]

    @property
    def ht_unten(self):
        if len(self.hts) == 2:
            return self.hts[0]
        else:
            return None

    def assert_checking(self):
        assert 1 <= len(self.hts) <= 2

        assert isinstance(self.ht_oben, Mantel)  # oben ist ein Mantel
        assert isinstance(self.kegel, Mantel)  # kegel ist ein Mantel
        assert self.kegel.d[0] < self.kegel.d[1]  # kegel ist tatsächlich ein mantel

        # h_oben, h_unten should really exist
        assert self.h_oben > 0
        assert self.h_unten > 0

        # if there is a hauptteil connecting from under
        if self.ht_unten is not None:
            assert len(self.hts) == 2
            assert isinstance(self.ht_unten, (Mantel, StandZarge))
            # assert self.r in [x * self.ht_oben.d[0] for x in [0.1, 0.154]]  # to make sure it can be calculated as STZ

    @property
    def pos_v(self):
        return self.ht_oben.pos_v

    @property
    def pos_h(self):
        return self.ht_oben.pos_h

    @classmethod
    def from_dict(cls, adict):
        try:
            kegelkrempe = KegelKrempe(name=adict['name'], h_non_kegelside=adict['h_oben'], h_kegelside=adict['h_unten'],
                                      r=adict['r'], hts=adict['hts'], kegel=adict['kegel'])
            assert all([hasattr(kegelkrempe, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating KegelKrempe: %s' % e)

        return kegelkrempe

    @property
    def make_segments(self):
        """
        this constructs the krempe
        constructs the arcs to that ht_oben and ht_unten are the tangents
        """

        if self.r > 0:
            # center of the arc
            _axis = LineSegment(Point((self.pos_h, self.pos_v)), Point((self.pos_h, self.pos_v + 1)))  # axis of the krempe
            # parallel to the upper segment
            l_oben = copy.deepcopy(self.ht_oben.segments[0])  # copy
            winkel_oben = _axis.angle_between_lines(l_oben)  # winkel between axis and segment of ht_oben
            l_oben.move(x0=-self.r / math.cos(math.radians(winkel_oben)), y0=0)  # moving to parallel
            # parallel to the lower segment
            l_unten = copy.deepcopy(self.kegel.segments[0])
            winkel_unten = _axis.angle_between_lines(l_unten)
            l_unten.move(x0=-self.r / math.cos(math.radians(winkel_unten)), y0=0)
            # center
            _mp = l_unten.intersection_point(l_oben)
            # the arc on the 'right' side. first point is on the lower ht, last is on the upper.
            _kre_re = arc(P=_mp, r=self.r, alpha_start=360-winkel_unten, alpha_end=360-winkel_oben, n=KEGELKREMPE_SEGMENTS)

            _graphrep_re = [Wall(copy.deepcopy(_kre_re.first_point), copy.deepcopy(_kre_re.last_point))]

        else:
            _kre_re = []
            _graphrep_re = []

        # calculating the new points at the ends of the abklinglange, in ht_oben and ht_unten
        _ls = [self.kegel.segments[0], self.ht_oben.segments[0]]
        _bbs = [self.kegel.bounding_box(), self.ht_oben.bounding_box()]

        # if h_oben and h_unten are to small, they are forced to yield points at the ends of the radii
        # however, if r=0, they only need to be positive
        _cut_v_pos = []
        if self.r > 0:
            # falls h_unten lang genug ist, wird das Kegel geschnitten
            if self.h_unten > abs(self.pos_v - _kre_re.first_point.y):
                _cut_v_pos.append(-max(self.h_unten, abs(self.pos_v - _kre_re.first_point.y)))
            else:
                _cut_v_pos.append(None)

            # falls h_oben lang genug ist, wird das Mantel über der Krempe geschnitten
            if self.h_oben > abs(self.pos_v - _kre_re.last_point.y):
                _cut_v_pos.append(max(self.h_oben, abs(self.pos_v - _kre_re.last_point.y)))
            else:
                _cut_v_pos.append(None)

        # falls r=0 wird geschnitten, falls die werte h_oben und h_unten nicht null sind
        elif self.r == 0:
            assert self.h_unten > 0
            assert self.h_oben > 0
            _cut_v_pos = [-self.h_unten, self.h_oben]

        else:
            raise Exception('WTF')

        _ips = []  # new points
        for line, bb, vpos in zip(_ls, _bbs, _cut_v_pos):
            if vpos is not None:
                _cutter = LineSegment(Point((min(bb[0])-abs(max(bb[1])), self.pos_v + vpos)),
                                      Point((max(bb[0])+abs(max(bb[1])), self.pos_v + vpos)))
                _ip = _cutter.intersection_point(line)
                _ips.append(_ip)
            else:
                _ips.append(None)

        # zugabe der segmente
        if self.r > 0:
            if _ips[1] is not None:
                _kre_re.append_point_at_end(_ips[1])
                _graphrep_re.append(Wall(_kre_re.segments[-1].i, _kre_re.segments[-1].j))

            if _ips[0] is not None:
                _kre_re._segments.insert(0, LineSegment(_ips[0], _kre_re.first_point))
                _graphrep_re.insert(0, Wall(_kre_re.segments[0].i, _kre_re.segments[0].j))

        else:
            _kre_re = PolygonLine([
                LineSegment(_ips[0], self.ht_oben.segments[0].i),
                LineSegment(self.kegel.segments[1].i, _ips[1])
            ])
            _graphrep_re = [Wall(*x.ij) for x in _kre_re.segments]

        _graphrep_li = copy.deepcopy(_graphrep_re)
        for w in _graphrep_li:
            w.mirror_vertical(x0=self.pos_h)
        self._graphrep = _graphrep_re + _graphrep_li

        # check: the endpoints of the arc MUST be on the lines of the upper and lower hts
        # checked is if this is true, in any order of ht_oben and ht_unten
        _fp = _kre_re.first_point
        _lp = _kre_re.last_point
        _e = (x.segments[0].is_point_on_line(_fp) and x.segments[0].is_internal_point(_fp) for x in [self.kegel, self.ht_oben])
        _z = (x.segments[0].is_point_on_line(_lp) and x.segments[0].is_internal_point(_lp) for x in [self.kegel, self.ht_oben])
        assert all(x or y for x, y in zip(_e, _z))

        # the arc on the 'left' side. This is done by copying the right side and mirroring
        _kre_li = copy.deepcopy(_kre_re)
        _kre_li.mirror_vertical(x0=self.pos_h)

        # modification of the ht above and the kegel
        # ht oben
        _ht = self.ht_oben
        p1 = Point(_kre_re.last_point.xy)
        p2 = Point((_ht.pos_h + _ht.d[1] / 2., _ht.pos_v + _ht.h))
        p3 = Point((_ht.pos_h - _ht.d[1] / 2., _ht.pos_v + _ht.h))
        p4 = Point(_kre_li.last_point.xy)
        ri_vert = Wall(p1, p2)
        le_vert = Wall(p3, p4)
        lo_hor = NonWall(p4, p1)
        hi_hor = NonWall(p2, p3)
        self.ht_oben._segments = [ri_vert, hi_hor, le_vert, lo_hor]
        self.ht_oben.graph.renew()

        _ht = self.kegel
        p1 = Point((_ht.pos_h + _ht.d[0] / 2., _ht.pos_v))
        p2 = Point(_kre_re.first_point.xy)
        p3 = Point(_kre_li.first_point.xy)
        p4 = Point((_ht.pos_h - _ht.d[0] / 2., _ht.pos_v))
        ri_vert = Wall(p1, p2)
        le_vert = Wall(p3, p4)
        lo_hor = NonWall(p4, p1)
        hi_hor = NonWall(p2, p3)
        self.kegel._segments = [ri_vert, hi_hor, le_vert, lo_hor]
        self.kegel.graph.renew()

        if self.ht_unten is not None:
            _ht = self.ht_unten
            p1 = Point((_ht.pos_h + _ht.d[0] / 2., _ht.pos_v))
            p2 = Point(_kre_re.last_segment.i.xy)
            p3 = Point(_kre_li.last_segment.i.xy)
            p4 = Point((_ht.pos_h - _ht.d[0] / 2., _ht.pos_v))
            ri_vert = Wall(p1, p2)
            le_vert = Wall(p3, p4)
            hi_hor = NonWall(p2, p3)
            if _ht.__class__.__name__ == 'Mantel':  # Mantel, Kegel
                lo_hor = NonWall(p4, p1)
            elif _ht.__class__.__name__ == 'StandZarge':
                lo_hor = Wall(p4, p1)
            else:
                raise Exception('WTF')
            _ht._segments = [ri_vert, hi_hor, le_vert, lo_hor]
            _ht.graph.renew()

        # self.ht_oben.plot(show=False)
        # self.kegel.plot(show=False)
        # self.ht_unten.plot(show=False)
        # _kre_re.plot(show=False)
        # _kre_li.plot(show=True)

        # preparing the return values
        if self.r > 0:
            _ret = [Wall(*x.ij) for x in _kre_re.segments] + \
                   [NonWall(_kre_re.last_point, _kre_li.last_point)] + \
                   [Wall(*x.ij) for x in _kre_li.segments] + \
                   [NonWall(_kre_li.first_point, _kre_re.first_point)]

        else:
            _ret = [Wall(*_kre_re.segments[0].ij), Wall(*_kre_re.segments[1].ij)] + \
                   [NonWall(_kre_re.segments[1].j, _kre_li.segments[1].j)] + \
                   [Wall(*_kre_li.segments[1].ij), Wall(*_kre_li.segments[0].ij)] + \
                   [NonWall(_kre_li.segments[0].i, _kre_re.segments[0].i)]

        # for r in [x for x in _ret if x.__class__.__name__ == 'Wall']:
        #     r.plot(show=False, style='k-')
        # for r in [x for x in _ret if x.__class__.__name__ == 'NonWall']:
        #     r.plot(show=False, style='k:')
        # r.plot(show=True, style='k:')

        return _ret

    @property
    def intex_points(self):
        """
        this is the SAME routine as the one for the class Mantel
        """

        _res = {}
        for s in self._graph.segments:
            _ret = []
            _res[s] = {'internal': [], 'external': [], 'indifferent': []}
            # the endpoints of a segment in ascending order vertically; lower, upper points
            _lp, _up = sorted(s.ij, key=lambda x: x.y)
            # "unit" for _lp, _up
            _e = [abs(min((_up.y - _lp.y) / 10., min((_p.x - self.pos_h) / 500., 20))) for _p in (_lp, _up)]
            # we create four points, two at each end of the segments
            for ratio in [0.01, 0.99]:
                _tp = point_on_line_at_ratio(line=s, ratio=ratio)
                _ret.append(Point((_tp.x + _e[0], _tp.y)))
                _ret.append(Point((_tp.x - _e[0], _tp.y)))

            for r in _ret:
                if self.is_internal_point(p=r) and not self.is_external_point(p=r):
                    _res[s]['internal'].append(r)
                elif not self.is_internal_point(p=r) and self.is_external_point(p=r):
                    _res[s]['external'].append(r)
                else:
                    raise Exception('Telling if Point internal or external impossible')

        #     self.plot(show=False, also_plot=_res[s]['internal'] + _res[s]['external'])
        # self.plot(show=True, also_plot=_res[s]['internal'] + _res[s]['external'])
        return _res

    @property
    def characteristic_size(self):
        return min(self.d) / 2.

    @property
    def isvalid(self):
        """
        A Hauptteil is valid if it is valid an its graph is valid.
        Since some subclasses are not valid polylines (e.g Flachboden) the evaluation is individual.
        """
        return True


class StandZarge(Mantel):
    # same as a Mantel but the lower segment is a wall
    def __init__(self, name=None, pos_v=None, pos_h=0, d=(), h=None):
        super(StandZarge, self).__init__(name=name, pos_v=pos_v, pos_h=pos_h, d=d, h=h)
        self._segments = self.make_segments

    @property
    def make_segments(self):
        _segs = super().make_segments
        _segs[3] = Wall(*_segs[3].ij)  # wird auf Wanll umgestellt
        return _segs


class RaumEnde(Hauptteil):
    def __init__(self, name=None, pos_v=None, d=None, pos_h=None):
        self.name = name
        self.d = d
        self.pos_h = pos_h
        self.pos_v = pos_v
        fb = Wall(Point((-self.d/2., 0)), Point((self.d/2., 0)))
        fb.move(x0=pos_h, y0=pos_v)
        super(RaumEnde, self).__init__(typ='raumende', segments=[fb])
        self._graph = Graph(hauptteil=self)
        self.graph.add_intex_points()

    @classmethod
    def from_dict(cls, adict):
        try:
            raumende = RaumEnde(
                name=adict['name'],
                pos_v=adict['pos_v'],
                pos_h=adict['pos_h'],
                d=adict['d'],
            )
            assert all([hasattr(raumende, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating Raumende: %s' % e)

        return raumende

    @property
    def intex_points(self):
        """
        RaumEnde hat keine intex punkte.
        
        Die Idee dahinter ist, dass diese bleiben unberechnet. Allerdings, soll es hier entschlossen werden?
        
        :return: 
        """
        _ret = {}
        for s in self.wall_segments:
            _ret[s] = {'internal': [], 'external': [], 'indifferent': []}
        return _ret

    @property
    def isvalid(self):
        """ valid if: valid as segment """
        return True


def run():

    re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
    kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
    ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
    bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
    mau = Mantel(name='mau', pos_h=0, pos_v=0, d=(4, 4), h=2)
    reu = RaumEnde(name='reu', pos_h=0, pos_v=0, d=4)
    kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0.4, hts=(ma, mau), kegel=kegel)
    beh19 = Behalter([ma, bo, re, kegel, kkr, mau, reu])
    beh19.plot(cycles=False)
    pp.pprint(beh19.raum_intex())

if __name__ == '__main__':

    print('')
    print('Anzahl Zeilen %d' % linenum())

    run()
