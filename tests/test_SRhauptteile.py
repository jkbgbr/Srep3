# -*- coding: utf-8 -*-

# from geometry.SR_line import LineSegment
# from geometry.SR_point import Point
# from geometry.SR_polygon import PolygonLine, square, convex_hull
from SR_Vessel.SR_hauptteile import Mantel, Korbbogenboden, Klopperboden, Flachboden, \
    Behalter, RaumEnde, Rohrboden, KegelKrempe, StandZarge, raum_dict_checker
import unittest
import math
from geometry import PURGE
import pprint as pp

EPS = 1e-10  # accepted numerical error


class TestPolygon_basics(unittest.TestCase):

    def setUp(self):
        # single tests
        self.ma1 = Mantel(pos_v=3, pos_h=4, d=(3, 6), h=2)
        self.fb1 = Flachboden(pos_v=3, pos_h=4, d=2)
        self.kb1 = Korbbogenboden(pos_v=3., pos_h=4., d=5., h1=0.1, orient=1)
        self.kb2 = Korbbogenboden(pos_v=3, pos_h=4, d=5, h1=0.1, orient=-1)
        self.kb1 = Korbbogenboden(pos_v=3000., pos_h=0., d=2500., h1=50, orient=1)
        self.ma11 = Mantel(pos_v=1000., pos_h=200., d=(1500, 2500), h=2000)
        self.fb11 = Flachboden(pos_v=500, pos_h=200, d=1350)

    def test_mantel_1(self):
        self.assertEqual(self.ma1.isvalid, True)

    def test_mantel_2(self):
        self.assertEqual(self.ma1.polygonline.isclosed, True)

    def test_mantel_3(self):
        self.assertEqual(self.ma1.polygonline.isconvex, True)

    def test_flachboden(self):
        self.assertEqual(self.fb1.isvalid, True)

    def test_korbbogenboden_1(self):
        self.assertEqual(self.kb1.isvalid, True)

    def test_korbbogenboden_2(self):
        self.assertEqual(self.kb2.isvalid, True)

    def test_internal_points_TorisphericalBoden_1(self):
        self.assertEqual(len(self.kb1.indifferent) == 0, True)

    def test_internal_points_TorisphericalBoden_2(self):
        self.assertEqual(all([self.kb1.graph.is_internal_point(p=P) for P in self.kb1.graph.internal]), True)

    def test_external_points_TorisphericalBoden_3(self):
        self.assertEqual(all([self.kb1.graph.is_external_point(p=P) for P in self.kb1.graph.external]), True)

    def test_external_points_Mantel_1(self):
        self.assertEqual(all([self.ma11.is_external_point(p=P) for P in self.ma11.graph.external]), True)

    def test_internal_points_Mantel_2(self):
        self.assertEqual(len(self.ma11.indifferent) == 0, True)

    def test_internal_points_Mantel_3(self):
        self.assertEqual(all([self.ma11.is_internal_point(p=P) for P in self.ma11.graph.internal]), True)

    def test_internal_points_Flachboden_1(self):
        self.assertEqual(len(self.fb11.internal) == 0, True)  # no internal points

    def test_external_points_Flachboden_2(self):
        self.assertEqual(len(self.fb11.external) == 0, True)  # no external points

    def test_internal_points_Flachboden_3(self):
        self.assertEqual(all([self.fb11.graph.is_internal_point(p=P) for P in self.fb11.graph.indifferent]), False)  # all indifferents are externals


class TestHauptteilCreation(unittest.TestCase):

    def tearDown(self):
        self.to = None
        self.adict = {}
        self.assertEqual(all([hasattr(self.to, x) for x in self.adict.keys()]), True)

    def test_mantel_from_dict(self):
        adict = {'name': 'mantel1', 'pos_v': 2, 'pos_h': 5, 'd': (3, 4), 'h': 3}
        self.to = Mantel.from_dict(adict)
        self.assertEqual(self.to.name, 'mantel1')
        self.assertEqual(self.to.pos_v, 2)
        self.assertEqual(self.to.pos_h, 5)
        self.assertEqual(self.to.d, (3, 4))
        self.assertEqual(self.to.h, 3)

    def test_standzarge_from_dict(self):
        adict = {'name': 'stz1', 'pos_v': 2, 'pos_h': 5, 'd': (3, 4), 'h': 3}
        self.to = Mantel.from_dict(adict)
        self.assertEqual(self.to.name, 'stz1')
        self.assertEqual(self.to.pos_v, 2)
        self.assertEqual(self.to.pos_h, 5)
        self.assertEqual(self.to.d, (3, 4))
        self.assertEqual(self.to.h, 3)

    def test_klopperboden_from_dict(self):
        adict = {'name': 'klopperboden1', 'pos_v': 2, 'pos_h': 5, 'd': 1, 'h1': 0.2, 'orient': 1}
        self.to = Klopperboden.from_dict(adict)
        self.assertEqual(self.to.name, 'klopperboden1')
        self.assertEqual(self.to.pos_v, 2)
        self.assertEqual(self.to.pos_h, 5)
        self.assertEqual(self.to.d, 1)
        self.assertEqual(self.to.h1, 0.2)
        self.assertEqual(self.to.orient, 1)

    def test_korbbogenboden_from_dict(self):
        adict = {'name': 'korbbogenboden1', 'pos_v': 2, 'pos_h': 5, 'd': 1, 'h1': 0.2, 'orient': -1}
        self.to = Korbbogenboden.from_dict(adict)
        self.assertEqual(self.to.name, 'korbbogenboden1')
        self.assertEqual(self.to.pos_v, 2)
        self.assertEqual(self.to.pos_h, 5)
        self.assertEqual(self.to.d, 1)
        self.assertEqual(self.to.h1, 0.2)
        self.assertEqual(self.to.orient, -1)

    def test_flachboden_from_dict(self):
        adict = {'name': 'flachboden1', 'pos_v': 2, 'pos_h': 5, 'd': 1}
        self.to = Flachboden.from_dict(adict)
        self.assertEqual(self.to.name, 'flachboden1')
        self.assertEqual(self.to.pos_v, 2)
        self.assertEqual(self.to.pos_h, 5)
        self.assertEqual(self.to.d, 1)

    # def test_rohrboden_from_dict(self):
    #     adict = {'name': 'rohrboden1', 'pos_v': 2, 'pos_h': 5, 'd': 1}
    #     rohrboden = Rohrboden.from_dict(adict)
    #     self.assertEqual(rohrboden.name, 'rohrboden1')
    #     self.assertEqual(rohrboden.pos_v, 2)
    #     self.assertEqual(rohrboden.pos_h, 5)
    #     self.assertEqual(rohrboden.d, 1)

    def test_kegelkrempe_from_dict(self):
        kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        ma = Mantel(name='ma_unten', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        adict = {'name': 'kegelkrempe1', 'h_oben': 0.4, 'h_unten': 0.3, 'r': 0.5, 'hts': (ma, kegel), 'kegel': kegel}
        Kegelkrempe = KegelKrempe.from_dict(adict)
        self.assertEqual(Kegelkrempe.name, 'kegelkrempe1')
        self.assertEqual(Kegelkrempe.r, 0.5)
        self.assertEqual(Kegelkrempe.h_oben, 0.4)
        self.assertEqual(Kegelkrempe.h_unten, 0.3)
        self.assertEqual(Kegelkrempe.ht_oben, ma)
        self.assertEqual(Kegelkrempe.ht_unten, kegel)


class TestPolygon_Behalter1(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # simple Behälter 1 test
        cls.ma22 = Mantel(name='ma22', pos_h=1, pos_v=0, d=(3, 3), h=2)
        cls.bo22 = Klopperboden(name='bo22', pos_h=1, pos_v=2, d=3, h1=0.1)
        cls.bu22 = Klopperboden(name='bu22', pos_h=1, pos_v=0, d=3, h1=0.1, orient=-1)
        cls.beh1 = Behalter(sorted([cls.ma22, cls.bo22, cls.bu22], key=lambda x: x.pos_v))

    def test_simple_behalter1_1(self):
        self.assertEqual(len(self.beh1.cycles.keys()) == 1, True)  # just one cycle

    def test_simple_behalter1_2(self):
        self.assertDictEqual(self.beh1.disjunct_cycles, {0: []}, True)

    def test_simple_behalter1_3(self):
        self.assertDictEqual(self.beh1.overlaps_in_cycles, {}, True)

    def test_simple_behalter1_4(self):
        self.assertDictEqual(self.beh1.contained_cycles, {0: []}, True)

    def test_simple_behalter1_5(self):
        self.assertTupleEqual(tuple([x.name for x in self.beh1.raume]), ('R1',), True)

    def test_simple_behalter1_6(self):
        for ht in self.beh1.polygons:
            self.assertDictEqual(self.beh1.raum_intex()[ht.name],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})


class TestPolygon_Behalter1_cut(unittest.TestCase):

    # todo: speed up the split algorithm
    @classmethod
    def setUpClass(cls):
        # simple Behälter 1 test. testing splits: the result must be the same
        cls.ma22 = Mantel(name='ma22', pos_h=1, pos_v=0, d=(3, 3), h=2)
        cls.bo22 = Klopperboden(name='bo22', pos_h=1, pos_v=2, d=3, h1=0.1)
        cls.bu22 = Klopperboden(name='bu22', pos_h=1, pos_v=0, d=3, h1=0.1, orient=-1)
        cls.ma22.split(v=-10.)
        cls.ma22.split(v=-0.1)
        cls.ma22.split(v=0.2)
        cls.ma22.split(v=1.0)
        cls.ma22.split(v=1.2)
        cls.ma22.split(v=1.9)
        cls.bo22.split(v=2.005)
        cls.bu22.split(v=-0.07)
        cls.bu22.split(v=-0.005)
        cls.beh1 = Behalter(sorted([cls.ma22, cls.bo22, cls.bu22], key=lambda x: x.pos_v))

    def test_split_fails(self):
        # Cut in Boden
        self.assertEqual(self.bo22.split(v=2.05), False)
        # # cut so that an endpoint is cut
        self.ma22.split(v=2)
        self.assertEqual(self.ma22.split(v=2), False)
        # # no intersection
        self.assertEqual(self.bo22.split(v=100), False)
        # # raumEnde is not cut
        self.re = RaumEnde(name='re', pos_h=1, pos_v=0, d=3)
        self.assertEqual(self.re.split(v=0), False)

    def test_simple_behalter1_1(self):
        self.assertEqual(len(self.beh1.graph.cycles.keys()) == 1, True)  # just one cycle
        self.assertEqual(len(self.beh1.cycles.keys()) == 1, True)  # just one cycle

    def test_simple_behalter1_2(self):
        self.assertDictEqual(self.beh1.disjunct_cycles, {0: []}, True)

    def test_simple_behalter1_3(self):
        self.assertDictEqual(self.beh1.overlaps_in_cycles, {}, True)

    def test_simple_behalter1_4(self):
        self.assertDictEqual(self.beh1.contained_cycles, {0: []}, True)

    def test_simple_behalter1_5(self):
        self.assertTupleEqual(tuple([x.name for x in self.beh1.raume]), ('R1',), True)

    def test_simple_behalter1_6(self):
        for ht in self.beh1.polygons:
            self.assertDictEqual(self.beh1.raum_intex()[ht.name],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})


class TestPolygon_Behalter2(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # simple Behälter 2 test
        cls.ma33 = Mantel(name='ma33', pos_h=1, pos_v=0, d=(3, 3), h=2)
        cls.fbo33 = Flachboden(name='fbo33', pos_h=1, pos_v=2, d=3)
        cls.fbu33 = Flachboden(name='fbu33', pos_h=1, pos_v=0, d=3)
        cls.beh2 = Behalter(sorted([cls.ma33, cls.fbo33, cls.fbu33], key=lambda x: x.pos_v))

    def test_simple_behalter2_1(self):
        self.assertEqual(len(self.beh2.cycles.keys()) == 1, True)  # just one cycle

    def test_simple_behalter2_2(self):
        self.assertDictEqual(self.beh2.disjunct_cycles, {0: []}, True)

    def test_simple_behalter2_3(self):
        self.assertDictEqual(self.beh2.overlaps_in_cycles, {}, True)

    def test_simple_behalter2_4(self):
        self.assertDictEqual(self.beh2.contained_cycles, {0: []}, True)

    def test_simple_behalter2_5(self):
        self.assertTupleEqual(tuple([x.name for x in self.beh2.raume]), ('R1',), True)

    def test_simple_behalter2_6(self):
        self.assertDictEqual(self.beh2.raum_intex()['ma33'],
                             {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})

    def test_simple_behalter2_7(self):
        for htname in ['fbo33', 'fbu33']:
            self.assertDictEqual(self.beh2.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': {'R1'}, 'innenraum': set()})


class TestPolygon_Behalter3(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # doppel Behälter 1 test: two raums next to each other
        cls.ma144 = Mantel(name='ma144', pos_h=1, pos_v=1, d=(3, 3), h=2)
        cls.ma244 = Mantel(name='ma244', pos_h=1, pos_v=3, d=(3, 3), h=2)
        cls.bo44 = Klopperboden(name='bo44', pos_h=1, pos_v=5, d=3, h1=0.1)
        cls.bm44 = Klopperboden(name='bm44', pos_h=1, pos_v=3, d=3, h1=0.1, orient=-1)
        cls.bu44 = Klopperboden(name='bu44', pos_h=1, pos_v=1, d=3, h1=0.1, orient=-1)
        cls.beh3 = Behalter(sorted([cls.ma144, cls.ma244, cls.bo44, cls.bm44, cls.bu44], key=lambda x: x.pos_v))

    def test_doppelraum_behalter_1_1(self):
        self.assertEqual(len(self.beh3.cycles.keys()) == 2, True)  # just one cycle

    def test_doppelraum_behalter_1_2(self):
        self.assertDictEqual(self.beh3.disjunct_cycles, {0: [], 1: []}, True)

    def test_doppelraum_behalter_1_3(self):
        self.assertDictEqual(self.beh3.overlaps_in_cycles_named, {(0, 1): {'bm44'}}, True)

    def test_doppelraum_behalter_1_4(self):
        self.assertDictEqual(self.beh3.contained_cycles, {0: [], 1: []}, True)

    def test_doppelraum_behalter_1_5(self):
        self.assertEqual(set([x.name for x in self.beh3.raume]), {'R1', 'R2'}, True)

    def test_doppelraum_behalter_1_6(self):
        _base = {'bm44': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'bo44': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'bu44': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma144': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma244': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(dict_to_comapare=self.beh3.raum_intex(), basedict=_base), True)


class TestPolygon_Behalter4(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # doppel Behälter 2 test: three raums next to each other
        cls.ma121 = Mantel(name='ma121', pos_v=0, pos_h=0, d=(4, 4), h=1)
        cls.ma122 = Mantel(name='ma122', pos_v=1, pos_h=0, d=(4, 4), h=1)
        cls.ma123 = Mantel(name='ma123', pos_v=2, pos_h=0, d=(4, 4), h=1)
        cls.bo12 = Korbbogenboden(name='bo12', pos_h=0, pos_v=3, d=4, h1=0.1)
        cls.bu12 = Klopperboden(name='bu12', pos_h=0, pos_v=0, d=4, h1=0.1, orient=-1)
        cls.fb121 = Flachboden(name='fb121', pos_v=1, d=4, pos_h=0)
        cls.fb122 = Flachboden(name='fb122', pos_v=2, d=4, pos_h=0)
        cls.beh4 = Behalter([cls.ma121, cls.ma122, cls.ma123, cls.bo12, cls.bu12, cls.fb121, cls.fb122])

    def test_doppelraum_behalter_2_1(self):
        self.assertEqual(len(self.beh4.cycles.keys()) == 3, True)

    def test_doppelraum_behalter_2_2(self):
        self.assertDictEqual(self.beh4.disjunct_cycles, {0: [], 1: [], 2: []}, True)

    # def test_doppelraum_behalter_2_3(self):
    #     self.assertDictEqual(self.beh4.overlaps_in_cycles_named, {(0, 1): {'fb121'}, (0, 2): {'fb122'}}, True)

    def test_doppelraum_behalter_2_4(self):
        self.assertDictEqual(self.beh4.contained_cycles, {0: [], 1: [], 2: []}, True)

    def test_doppelraum_behalter_2_5(self):
        self.assertEqual(set([x.name for x in self.beh4.raume]), {'R1', 'R2', 'R3'}, True)

    def test_doppelraum_behalter_2_6(self):
        _base = {'bo12': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bu12': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R3'}},
                 'fb121': {'aussenraum': set(), 'indifferent': {'R1', 'R3'}, 'innenraum': set()},
                 'fb122': {'aussenraum': set(), 'indifferent': {'R1', 'R2'}, 'innenraum': set()},
                 'ma121': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R3'}},
                 'ma122': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma123': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh4.raum_intex()), True)


class TestPolygon_Behalter5(unittest.TestCase):
    """ Doppelmantelbehälter """

    @classmethod
    def setUpClass(cls):
        cls.ma341 = Mantel(name='ma341', pos_v=0, pos_h=0, d=(4000, 4000), h=1000)
        cls.bo343 = Klopperboden(name='bo343', pos_h=0, pos_v=1000, d=4000, h1=100)
        cls.bu343 = Klopperboden(name='bu343', pos_h=0, pos_v=0, d=4000, h1=100, orient=-1)
        cls.ma341a = Mantel(name='ma341a', pos_v=-100, pos_h=0, d=(4100, 4100), h=1200)
        cls.bo343a = Klopperboden(name='bo343a', pos_h=0, pos_v=1100, d=4100, h1=100)
        cls.bu343a = Klopperboden(name='bu343a', pos_h=0, pos_v=-100, d=4100, h1=100, orient=-1)
        cls.beh5 = Behalter(sorted([cls.ma341, cls.bo343, cls.bu343, cls.ma341a, cls.bo343a, cls.bu343a], key=lambda x: x.pos_v))
        # cls.beh5.plot(cycles=True)

    def test_doppelmantel_behalter_1(self):
        self.assertEqual(len(self.beh5.cycles.keys()) == 2, True)

    def test_doppelmantel_behalter_1_2(self):
        self.assertDictEqual(self.beh5.disjunct_cycles, {0: [], 1: []}, True)

    def test_doppelmantel_behalter_1_3(self):
        self.assertDictEqual(self.beh5.overlaps_in_cycles_named, {}, True)

    def test_doppelmantel_behalter_1_4(self):
        self.assertDictEqual(self.beh5.contained_cycles, {0: [1], 1: []}, True)

    def test_doppelmantel_behalter_1_5(self):
        self.assertEqual(set([x.name for x in self.beh5.raume]), {'R1', 'R2'}, True)

    def test_doppelmantel_behalter_1_7(self):
        _base = {'bo343': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'bo343a': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'bu343': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'bu343a': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma341': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma341a': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh5.raum_intex()), True)


class TestPolygon_Behalter5_cut(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma341 = Mantel(name='ma341', pos_v=0, pos_h=0, d=(4000, 4000), h=1000)
        cls.bo343 = Klopperboden(name='bo343', pos_h=0, pos_v=1000, d=4000, h1=100)
        cls.bu343 = Klopperboden(name='bu343', pos_h=0, pos_v=0, d=4000, h1=100, orient=-1)
        cls.ma341a = Mantel(name='ma341a', pos_v=-100, pos_h=0, d=(4100, 4100), h=1200)
        cls.bo343a = Klopperboden(name='bo343a', pos_h=0, pos_v=1100, d=4100, h1=100)
        cls.bu343a = Klopperboden(name='bu343a', pos_h=0, pos_v=-100, d=4100, h1=100, orient=-1)
        cls.ma341.split(v=-2000)
        cls.ma341.split(v=800)
        cls.ma341a.split(v=700)
        cls.beh5 = Behalter(sorted([cls.ma341, cls.bo343, cls.bu343, cls.ma341a, cls.bo343a, cls.bu343a], key=lambda x: x.pos_v))

    def test_doppelmantel_behalter_1(self):
        self.assertEqual(len(self.beh5.cycles.keys()) == 2, True)  # just one cycle

    def test_doppelmantel_behalter_1_2(self):
        self.assertDictEqual(self.beh5.disjunct_cycles, {0: [], 1: []}, True)

    def test_doppelmantel_behalter_1_3(self):
        self.assertDictEqual(self.beh5.overlaps_in_cycles_named, {}, True)

    def test_doppelmantel_behalter_1_4(self):
        self.assertDictEqual(self.beh5.contained_cycles, {0: [1], 1: []}, True)

    def test_doppelmantel_behalter_1_5(self):
        self.assertEqual(set([x.name for x in self.beh5.raume]), {'R1', 'R2'}, True)

    def test_doppelmantel_behalter_1_6(self):
        for htname in ['bo343', 'bu343', 'ma341']:
            self.assertDictEqual(self.beh5.raum_intex()[htname],
                                 {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}})

    def test_doppelmantel_behalter_1_7(self):
        for htname in ['bo343a', 'bu343a', 'ma341a']:
            self.assertDictEqual(self.beh5.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})


class TestPolygon_Behalter6(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.ma61 = Mantel(name='ma61', pos_v=-2, pos_h=5, d=(4, 4), h=10)
        cls.bu61a = Klopperboden(name='bu61a', pos_v=-2, pos_h=5, d=4, h1=0.1, orient=-1)
        cls.bo61a = Klopperboden(name='bo61a', pos_v=8, pos_h=5, d=4, h1=0.1)

        cls.bu62a = Klopperboden(name='bu62a', pos_v=-1, pos_h=5, d=3.9, h1=0.1, orient=-1)
        cls.ma61a = Mantel(name='ma61a', pos_v=-1.0, pos_h=5, d=(3.9, 3.9), h=1)
        cls.fb62b = Flachboden(name='fb62b', pos_v=0, pos_h=5, d=3.9)

        cls.fb62c = Flachboden(name='fb62c', pos_v=1, pos_h=5, d=3.9)
        cls.ma61b = Mantel(name='ma61b', pos_v=1, pos_h=5, d=(3.9, 3.9), h=1)
        cls.fb62d = Flachboden(name='fb62d', pos_v=2, pos_h=5, d=3.9)

        cls.ma61e = Mantel(name='ma61e', pos_v=0, pos_h=5, d=(3.9, 3.9), h=1)
        cls.ma63f = Mantel(name='ma63f', pos_v=2, pos_h=5, d=(3.9, 3.9), h=0.1)
        cls.bo63e = Klopperboden(name='bo63e', pos_v=2.1, pos_h=5, d=3.9, h1=1)

        cls.bu163d = Klopperboden(name='bu163d', pos_v=5, pos_h=5, d=3.9, h1=0.1, orient=-1)
        cls.ma161d = Mantel(name='ma161d', pos_v=5, pos_h=5, d=(3.9, 3.9), h=1)
        cls.bo163e = Klopperboden(name='bo163e', pos_v=6, pos_h=5, d=3.9, h1=0.1)

        cls.beh6 = Behalter(sorted([cls.ma61, cls.ma61a, cls.ma61b, cls.bu61a, cls.bo61a, cls.bu62a,
                                    cls.fb62b, cls.fb62d, cls.fb62c, cls.ma61e, cls.ma63f, cls.bo63e, cls.bu163d,
                                    cls.ma161d, cls.bo163e], key=lambda x: x.pos_v))
        # cls.beh6.plot(cycles=False)
#
    def test_complexmantel_behalter_1(self):
        self.assertEqual(len(self.beh6.cycles.keys()) == 6, True)

    def test_complexmantel_behalter_1_2(self):
        self.assertDictEqual(self.beh6.disjunct_cycles, {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}, True)

    # def test_complexmantel_behalter_1_3(self):
    #     self.assertDictEqual(self.beh6.overlaps_in_cycles_named, {(0, 1): {'fb62c'}, (1, 3): {'fb62d'}, (0, 2): {'fb62b'}}, True)
    #
    # def test_complexmantel_behalter_1_4(self):
    #     self.assertDictEqual(self.beh6.contained_cycles, {0: [], 1: [], 2: [], 3: [], 4: [0, 1, 2, 3, 5], 5: []}, True)

    def test_complexmantel_behalter_1_5(self):
        self.assertEqual(set([x.name for x in self.beh6.raume]), {'R1', 'R2', 'R3', 'R4', 'R5', 'R6'}, True)

    def test_complexmantel_behalter_1_6(self):
        _base = {'bo163e': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R5'}},
                 'bo61a': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R4'}},
                 'bo63e': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R3'}},
                 'bu163d': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R5'}},
                 'bu61a': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R4'}},
                 'bu62a': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R6'}},
                 'fb62b': {'aussenraum': set(), 'indifferent': {'R6', 'R2'}, 'innenraum': set()},
                 'fb62c': {'aussenraum': set(), 'indifferent': {'R2', 'R1'}, 'innenraum': set()},
                 'fb62d': {'aussenraum': set(), 'indifferent': {'R3', 'R1'}, 'innenraum': set()},
                 'ma161d': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R5'}},
                 'ma61': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R4'}},
                 'ma61a': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R6'}},
                 'ma61b': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma61e': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma63f': {'aussenraum': {'R4'}, 'indifferent': set(), 'innenraum': {'R3'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh6.raum_intex()), True)


class TestPolygon_Behalter7(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma211 = Mantel(name='mau', pos_v=0, pos_h=0, d=(4, 4), h=2)
        cls.ma212 = Mantel(name='mao', pos_v=2, pos_h=0, d=(4, 4), h=2)
        cls.bo221 = Klopperboden(name='bo', pos_h=0, pos_v=4, d=4, h1=0.1)
        cls.bu221 = Klopperboden(name='bu', pos_h=0, pos_v=0, d=4, h1=0.1, orient=-1)
        cls.zb221 = Klopperboden(name='zb', pos_h=0, pos_v=2, d=4, h1=0.1, orient=-1)
        cls.beh7 = Behalter(sorted([cls.ma211, cls.ma212, cls.bo221, cls.bu221, cls.zb221], key=lambda x: x.pos_v))

    def test_zweiraummantel_behalter_7(self):
        self.assertEqual(len(self.beh7.cycles.keys()) == 2, True)

    def test_zweiraummantel_behalter_7_2(self):
        self.assertDictEqual(self.beh7.disjunct_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_7_3(self):
        self.assertDictEqual(self.beh7.overlaps_in_cycles_named, {(0, 1): {'zb'}}, True)

    def test_zweiraummantel_behalter_7_4(self):
        self.assertDictEqual(self.beh7.contained_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_7_5(self):
        self.assertEqual(set([x.name for x in self.beh7.raume]), {'R1', 'R2'}, True)

    def test_zweiraummantel_behalter_7_6(self):
        _base = {'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'bu': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'mao': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'mau': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'zb': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh7.raum_intex()), True)

#
#
# class TestPolygon_Behalter7_cut(unittest.TestCase):
#
    # @classmethod
#     def setUpClass(cls):
#         cls.ma211 = Mantel(name='mau', pos_v=0, pos_h=0, d=(4, 4), h=2)
#         cls.ma212 = Mantel(name='mao', pos_v=2, pos_h=0, d=(4, 4), h=2)
#         cls.bo221 = Klopperboden(name='bo', pos_h=0, pos_v=4, d=4, h1=0.1)
#         cls.bu221 = Klopperboden(name='bu', pos_h=0, pos_v=0, d=4, h1=0.1, orient=-1)
#         cls.zb221 = Klopperboden(name='zb', pos_h=0, pos_v=2, d=4, h1=0.1, orient=-1)
#         cls.ma211.split(v=1.2)
#         cls.ma212.split(v=3.1)
#         cls.bo221.split(v=4.06)
#         cls.bu221.split(v=-0.04)
#         cls.zb221.split(v=2-0.04)
#         cls.beh7 = Behalter(sorted([cls.ma211, cls.ma212, cls.bo221, cls.bu221, cls.zb221], key=lambda x: x.pos_v))
#
#     def test_zweiraummantel_behalter_7(self):
#         self.assertEqual(len(self.beh7.cycles.keys()) == 2, True)
#
#     def test_zweiraummantel_behalter_7_2(self):
#         self.assertDictEqual(self.beh7.disjunct_cycles, {0: [], 1: []}, True)
#
#     def test_zweiraummantel_behalter_7_3(self):
#         self.assertDictEqual(self.beh7.overlaps_in_cycles_named, {(0, 1): {'zb'}}, True)
#
#     def test_zweiraummantel_behalter_7_4(self):
#         self.assertDictEqual(self.beh7.contained_cycles, {0: [], 1: []}, True)
#
#     def test_zweiraummantel_behalter_7_5(self):
#         self.assertEqual(set([x.name for x in self.beh7.raume]), {'R1', 'R2'}, True)
#
#     def test_zweiraummantel_behalter_7_6(self):
#         for htname in ['mao', 'bo']:
#             self.assertDictEqual(self.beh7.raum_intex()[htname],
#                                  {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}})
#
#     def test_zweiraummantel_behalter_7_7(self):
#         for htname in ['mau', 'bu']:
#             self.assertDictEqual(self.beh7.raum_intex()[htname],
#                                  {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})
#
#     def test_zweiraummantel_behalter_7_8(self):
#         self.assertDictEqual(self.beh7.raum_intex()['zb'], {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}})
#
#
class TestPolygon_Behalter8(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma331 = Mantel(name='ma331', pos_v=1.5, pos_h=0, d=(4, 4), h=1.7)
        cls.ma332 = Mantel(name='ma332', pos_v=1, pos_h=0, d=(4.2, 4.2), h=1.2)
        cls.bui33 = Korbbogenboden(name='bui33', pos_h=0, pos_v=1.5, d=4, h1=0.1, orient=-1)
        cls.bua33 = Klopperboden(name='bua33', pos_h=0, pos_v=1, d=4.2, h1=0.1, orient=-1)
        cls.re331 = RaumEnde(name='re331', pos_v=2.2, d=4.2, pos_h=0)
        cls.fb332 = Klopperboden(name='fb332', pos_v=3.2, d=4, pos_h=0, h1=0.1)
        cls.beh8 = Behalter([cls.ma331, cls.ma332, cls.bui33, cls.bua33, cls.re331, cls.fb332])
        # cls.beh8.plot(cycles=True)

    def test_zweiraummantel_behalter_8(self):
        self.assertEqual(len(self.beh8.cycles.keys()) == 2, True)

    def test_zweiraummantel_behalter_8_2(self):
        self.assertDictEqual(self.beh8.disjunct_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_8_3(self):
        self.assertDictEqual(self.beh8.overlaps_in_cycles_named, {(0, 1): {'ma331', 'bui33'}}, True)

    def test_zweiraummantel_behalter_8_4(self):
        self.assertDictEqual(self.beh8.contained_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_8_5(self):
        self.assertEqual(set([x.name for x in self.beh8.raume]), {'R1', 'R2'}, True)

    def test_zweiraummantel_behalter_8_6(self):
        for htname in ['bua33', 'ma332']:
            self.assertDictEqual(self.beh8.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}})

    def test_zweiraummantel_behalter_8_7(self):
        for htname in ['fb332']:
            self.assertDictEqual(self.beh8.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_8_8(self):
        for htname in ['bui33', 'ma331']:
            self.assertDictEqual(self.beh8.raum_intex()[htname],
                                 {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_8_9(self):
        _base = {'bua33': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bui33': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'fb332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma331': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're331': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh8.raum_intex()), True)


class TestPolygon_Behalter9(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma331 = Mantel(name='ma331', pos_v=1.5, pos_h=0, d=(4, 4), h=1.7)
        cls.ma332 = Mantel(name='ma332', pos_v=1, pos_h=0, d=(4.2, 4.2), h=1.2)
        cls.ma333 = Mantel(name='ma333', pos_v=3.2, pos_h=0, d=(4, 4), h=1.2)
        cls.bui33 = Korbbogenboden(name='bui33', pos_h=0, pos_v=1.5, d=4, h1=0.1, orient=-1)
        cls.bua33 = Klopperboden(name='bua33', pos_h=0, pos_v=1, d=4.2, h1=0.1, orient=-1)
        cls.re331 = RaumEnde(name='re331', pos_v=2.2, d=4.2, pos_h=0)
        cls.fb332 = Klopperboden(name='fb332', pos_v=4.4, d=4, pos_h=0, h1=0.1)
        cls.beh9 = Behalter([cls.ma331, cls.ma332, cls.ma333, cls.bui33, cls.bua33, cls.re331, cls.fb332])
        # cls.beh9.plot(cycles=True)

    def test_zweiraummantel_behalter_9(self):
        self.assertEqual(len(self.beh9.cycles.keys()) == 2, True)
#
    def test_zweiraummantel_behalter_9_2(self):
        self.assertDictEqual(self.beh9.disjunct_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_9_3(self):
        self.assertDictEqual(self.beh9.overlaps_in_cycles_named, {(0, 1): {'ma331', 'bui33'}}, True)

    def test_zweiraummantel_behalter_9_4(self):
        self.assertDictEqual(self.beh9.contained_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_9_5(self):
        self.assertEqual(set([x.name for x in self.beh9.raume]), {'R1', 'R2'}, True)

    def test_zweiraummantel_behalter_9_6(self):
        for htname in ['bua33', 'ma332']:
            self.assertDictEqual(self.beh9.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}})

    def test_zweiraummantel_behalter_9_7(self):
        for htname in ['fb332', 'ma333']:
            self.assertDictEqual(self.beh9.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_9_8(self):
        for htname in ['bui33', 'ma331']:
            self.assertDictEqual(self.beh9.raum_intex()[htname],
                                 {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_9_9(self):
        _base = {'bua33': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bui33': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'fb332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma333': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma331': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're331': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh9.raum_intex()), True)


class TestPolygon_Behalter10(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma331 = Mantel(name='ma331', pos_v=1.5, pos_h=0, d=(4, 4), h=1.7)
        cls.ma333 = Mantel(name='ma333', pos_v=3.2, pos_h=0, d=(4, 4), h=1.2)
        cls.ma332 = Mantel(name='ma332', pos_v=1, pos_h=0, d=(4.2, 4.2), h=2.3)
        cls.bui33 = Korbbogenboden(name='bui33', pos_h=0, pos_v=1.5, d=4, h1=0.1, orient=-1)
        cls.bua33 = Klopperboden(name='bua33', pos_h=0, pos_v=1, d=4.2, h1=0.1, orient=-1)
        cls.re331 = RaumEnde(name='re331', pos_v=3.3, d=4.2, pos_h=0)
        cls.fb332 = Klopperboden(name='fb332', pos_v=4.4, d=4, pos_h=0, h1=0.1)
        cls.beh10 = Behalter([cls.ma331, cls.ma332, cls.ma333, cls.bui33, cls.bua33, cls.re331, cls.fb332])
        # cls.beh10.plot(cycles=True)

    def test_zweiraummantel_behalter_9(self):
        self.assertEqual(len(self.beh10.cycles.keys()) == 2, True)

    def test_zweiraummantel_behalter_9_2(self):
        self.assertDictEqual(self.beh10.disjunct_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_9_3(self):
        print(self.beh10.overlaps_in_cycles_named)
        self.assertDictEqual(self.beh10.overlaps_in_cycles_named, {(0, 1): {'ma331', 'bui33', 'ma333'}}, True)

    def test_zweiraummantel_behalter_9_4(self):
        self.assertDictEqual(self.beh10.contained_cycles, {0: [], 1: []}, True)

    def test_zweiraummantel_behalter_9_5(self):
        self.assertEqual(set([x.name for x in self.beh10.raume]), {'R1', 'R2'}, True)

    def test_zweiraummantel_behalter_9_6(self):
        for htname in ['bua33', 'ma332']:
            self.assertDictEqual(self.beh10.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}})

    def test_zweiraummantel_behalter_9_7(self):
        for htname in ['fb332']:
            self.assertDictEqual(self.beh10.raum_intex()[htname],
                                 {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_9_8(self):
        for htname in ['bui33', 'ma331', 'ma333']:
            self.assertDictEqual(self.beh10.raum_intex()[htname],
                                 {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}})

    def test_zweiraummantel_behalter_9_9(self):
        _base = {'bua33': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bui33': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'fb332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma333': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma331': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma332': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're331': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh10.raum_intex()), True)


class TestPolygon_Behalter11(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ma_in = Mantel(name='ma_in', pos_h=0, pos_v=1, d=(4, 4), h=2.5)
        cls.bu_in = Korbbogenboden(name='bu_in', pos_h=0, pos_v=1, d=4, h1=0.1, orient=-1)
        cls.re_in = RaumEnde(name='re_in', pos_v=3.5, d=4, pos_h=0)
        cls.ma_au = Mantel(name='ma_au', pos_h=0, pos_v=0.9, d=(4.1, 4.1), h=2.1)
        cls.bu_au = Korbbogenboden(name='bu_au', pos_v=0.9, pos_h=0, d=4.1, h1=0.1, orient=-1)
        cls.re_au = RaumEnde(name='re_au', pos_v=3, d=4.1, pos_h=0)
        cls.beh11 = Behalter([cls.ma_in, cls.ma_au, cls.bu_in, cls.bu_au, cls.re_in, cls.re_au])
        # cls.beh11.plot(cycles=True)

    def test_zweiraummantel_behalter_1(self):
        self.assertEqual(len(self.beh11.cycles.keys()) == 2, True)

    def test_zweiraummantel_behalter_2(self):
        _base = {'bu_au': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bu_in': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma_au': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma_in': {'aussenraum': {'R2'}, 'indifferent': set(), 'innenraum': {'R1'}},
                 're_au': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 're_in': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh11.raum_intex()), True)

class TestPolygon_Behalter12(unittest.TestCase):

    # doppelstockbehalter mit STZ

    @classmethod
    def setUpClass(cls):
        # untere Teile
        cls.mau = Mantel(name='ma_unten', pos_h=0, pos_v=1, d=(4, 4), h=2.5)
        cls.mau.split(v=2)
        cls.buu = Korbbogenboden(name='bu_unten', pos_h=0, pos_v=1, d=4, h1=0.1, orient=-1)
        cls.buo = Korbbogenboden(name='bo_unten', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)

        # zarge
        cls.zarge = Mantel(name='zarge', pos_h=0, pos_v=3.5, d=(4, 4), h=2.5)

        # oben
        cls.mao = Mantel(name='ma_oben', pos_h=0, pos_v=6, d=(4, 4), h=2.5)
        cls.bou = Korbbogenboden(name='bu_oben', pos_h=0, pos_v=6, d=4, h1=0.1, orient=-1)
        cls.boo = Korbbogenboden(name='bo_oben', pos_h=0, pos_v=8.5, d=4, h1=0.1, orient=1)

        # standzarge
        cls.stz = StandZarge(name='stz', pos_h=0, pos_v=-1, d=(4, 4), h=2)

        cls.beh12 = Behalter([cls.mau, cls.buu, cls.buo, cls.mao, cls.bou, cls.boo, cls.zarge, cls.stz])
        # cls.beh12.plot(cycles=True)

    def test_doppelstockbehalter_behalter_1(self):
        self.assertEqual(len(self.beh12.cycles.keys()) == 4, True)

    def test_doppelstockbehalter_behalter_2(self):
        _base = {'bo_oben': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R4'}},
                 'bo_unten': {'aussenraum': {'R3'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'bu_oben': {'aussenraum': {'R3'}, 'indifferent': set(), 'innenraum': {'R4'}},
                 'bu_unten': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma_oben': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R4'}},
                 'ma_unten': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'stz': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'zarge': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R3'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh12.raum_intex()), True)


class TestPolygon_Behalter13(unittest.TestCase):

    # Behalter mit einem Raum und STZ. Viele schüsse im Mantel
    @classmethod
    def setUpClass(cls):
        # untere Teile
        cls.mau = Mantel(name='ma_unten', pos_h=0, pos_v=1, d=(4, 4), h=8.5)
        cls.buu = Korbbogenboden(name='bu_unten', pos_h=0, pos_v=1, d=4, h1=0.1, orient=-1)
        cls.buo = Korbbogenboden(name='bo_unten', pos_h=0, pos_v=9.5, d=4, h1=0.1, orient=1)
        cls.stz = StandZarge(name='stz', pos_h=0, pos_v=-1, d=(4, 4), h=2)
        cls.beh13 = Behalter([cls.mau, cls.buu, cls.buo, cls.stz])
        # cls.beh13.plot(cycles=True)

    def test_doppelstockbehalter_behalter_1(self):
        self.assertEqual(len(self.beh13.cycles.keys()) == 2, True)

    def test_doppelstockbehalter_behalter_2(self):
        _base = {'bo_unten': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'bu_unten': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma_unten': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'stz': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh13.raum_intex()), True)


class TestPolygon_Behalter14(unittest.TestCase):

    # Behalter mit kegel und kegelkrempe. Kegelkrempe hat r>0
    @classmethod
    def setUpClass(cls):
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0.5, hts=(cls.ma,), kegel=cls.kegel)
        cls.beh14 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr])
        # cls.beh15.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh14.cycles.keys()) == 1, True)

    def test_kegelbehalter_behalter_2(self):
        _base = {'re': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 'kegel': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'kkr': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh14.raum_intex()), True)


class TestPolygon_Behalter15(unittest.TestCase):

    # Behalter mit kegel und kegelkrempe. Krempe hat r=0
    @classmethod
    def setUpClass(cls):
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0., hts=(cls.ma,), kegel=cls.kegel)
        cls.beh15 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr])
        # cls.beh15.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh15.cycles.keys()) == 1, True)

    def test_kegelbehalter_behalter_2(self):
        _base = {'re': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 'kegel': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'kkr': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh15.raum_intex()), True)


class TestPolygon_Behalter16(unittest.TestCase):

    # Doppelmantelbehalter mit kegel und kegelkrempe. eine Krempe hat r=0, andere r>0
    @classmethod
    def setUpClass(cls):
        # outer
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0.5, hts=(cls.ma,), kegel=cls.kegel)

        # inner
        cls.rei = RaumEnde(name='re', pos_h=0, pos_v=1.1, d=1.8)
        cls.kegeli = Mantel(name='kegel', pos_h=0, pos_v=1.1, d=(1.8, 3.8), h=0.9)
        cls.mai = Mantel(name='ma', pos_h=0, pos_v=2, d=(3.8, 3.9), h=1.4)
        cls.boi = Korbbogenboden(name='bo', pos_h=0, pos_v=3.4, d=3.9, h1=0.1, orient=1)
        cls.kkri = KegelKrempe(name='kkr', h_non_kegelside=0.2, h_kegelside=0.15, r=0., hts=(cls.mai,), kegel=cls.kegeli)

        cls.beh16 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr,
                              cls.mai, cls.boi, cls.rei, cls.kegeli, cls.kkri])
        # cls.beh16.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh16.cycles.keys()) == 2, True)

    # def test_kegelbehalter_behalter_2(self):
    #     self.assertDictEqual(self.beh15.raum_intex(),
    #                          {'re': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
    #                           'kegel': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
    #                           'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
    #                           'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
    #                           'kkr': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
    #                          )

class TestPolygon_Behalter17(unittest.TestCase):

    # Behalter mit kegelkrempe uns standzarge. Krempe hat r=0
    @classmethod
    def setUpClass(cls):
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.stz = StandZarge(name='stz', pos_h=0, pos_v=0, d=(4, 4), h=2)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0., hts=(cls.ma, cls.stz), kegel=cls.kegel)
        cls.beh17 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr, cls.stz])
        # cls.beh19.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh17.cycles.keys()) == 2, True)

    def test_kegelbehalter_behalter_2(self):
        _base = {'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'kegel': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'kkr': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 'stz': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh17.raum_intex()), True)


class TestPolygon_Behalter18(unittest.TestCase):

    # Behalter mit kegelkrempe und standzarge. Krempe hat r>0
    @classmethod
    def setUpClass(cls):
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.stz = StandZarge(name='stz', pos_h=0, pos_v=0, d=(4, 4), h=2)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0.5, hts=(cls.ma, cls.stz), kegel=cls.kegel)
        cls.beh18 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr, cls.stz])
        # cls.beh19.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh18.cycles.keys()) == 2, True)

    def test_kegelbehalter_behalter_2(self):
        _base = {'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'kegel': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'kkr': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 'stz': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh18.raum_intex()), True)


class TestPolygon_Behalter19(unittest.TestCase):

    # Behalter mit kegelkrempe uns standzarge. Krempe hat r>0
    @classmethod
    def setUpClass(cls):
        cls.re = RaumEnde(name='re', pos_h=0, pos_v=1, d=2)
        cls.kegel = Mantel(name='kegel', pos_h=0, pos_v=1, d=(2, 4), h=1)
        cls.ma = Mantel(name='ma', pos_h=0, pos_v=2, d=(4, 4), h=1.5)
        cls.bo = Korbbogenboden(name='bo', pos_h=0, pos_v=3.5, d=4, h1=0.1, orient=1)
        cls.mau = Mantel(name='mau', pos_h=0, pos_v=0, d=(4, 4), h=2)
        cls.reu = RaumEnde(name='reu', pos_h=0, pos_v=0, d=4)
        cls.kkr = KegelKrempe(name='kkr', h_non_kegelside=0.4, h_kegelside=0.3, r=0.4, hts=(cls.ma, cls.mau), kegel=cls.kegel)
        cls.beh19 = Behalter([cls.ma, cls.bo, cls.re, cls.kegel, cls.kkr, cls.mau, cls.reu])
        # cls.beh19.plot(cycles=True)

    def test_kegelbehalter_behalter_1(self):
        self.assertEqual(len(self.beh19.cycles.keys()) == 2, True)

    def test_kegelbehalter_behalter_2(self):
        pp.pprint(self.beh19.raum_intex())
        _base = {'bo': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 'kegel': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'kkr': {'aussenraum': {'R1'}, 'indifferent': set(), 'innenraum': {'R2'}},
                 'ma': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R2'}},
                 're': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()},
                 'mau': {'aussenraum': set(), 'indifferent': set(), 'innenraum': {'R1'}},
                 'reu': {'aussenraum': set(), 'indifferent': set(), 'innenraum': set()}}
        self.assertEqual(raum_dict_checker(basedict=_base, dict_to_comapare=self.beh19.raum_intex()), True)
