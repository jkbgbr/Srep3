# -*- coding: utf-8 -*-

import math
from cached_property import cached_property
from geometry import _plotting_available, plt

EMPTY = (None, None)
PRECISION = 8  # that many decimals. rounding to this numer of decimals.
# needed to fix number representation issues when doing floating arithmetic
EPS = 1. / (10 ** PRECISION)  # accepted numerical error


def distance(p1=EMPTY, p2=EMPTY):
    """ The distance of two point_set """
    assert p1.isvalid
    assert p2.isvalid

    return math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)


class Point(object):
    def __init__(self, args=EMPTY):
        """
        creates a Point object in the 2D plane defined by a tuple
        """
        # if not args:
        #     args = EMPTY
        self._set_coords(coords=args)
        self.plot_style = None
        # self.ID = uuid.uuid4()

    def __repr__(self):
        return 'Point((%r, %r))' % (self.x, self.y)

    def __str__(self):
        if self.isvalid:
            return self.__repr__()
        else:
            return 'Non-valid instance of class %s' % self.__class__.__name__

    def __sub__(self, other):
        if self.isvalid:
            return Point((self.x - other.x, self.y - other.y))
        else:
            return 'Non-valid instance of class %s' % self.__class__.__name__

    def __add__(self, other):
        if self.isvalid:
            return Point((self.x + other.x, self.y + other.y))
        else:
            return 'Non-valid instance of class %s' % self.__class__.__name__

    def __eq__(self, other):
        if not self.isvalid or not other.isvalid:
            return False

        # if self.xy == other.xy and not self.isvalid:  # none are valid - still equal :-)
        #     return True
        # if not self.isvalid or not other.isvalid:
        #     return False

        return (
            type(other) == type(self) and
            all([abs(x-y) < EPS * 1e3 for x, y in zip(self.xy, other.xy)])
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self._coords)

    @property
    def x(self):
        """ Return x coordinate. """
        return self._coords[0]

    @property
    def y(self):
        """ Return y coordinate. """
        return self._coords[1]

    @property
    def xy(self):
        """Return x, y coordinates as a tuple."""
        return self.x, self.y

    def _set_coords(self, coords=EMPTY):
        """ sets coordinates; use it to overwrite current coordinates """
        try:
            # making sure floats are used, even if defined in int (Python 2)
            self._coords = tuple([round(float(x), PRECISION) for x in coords])
        except TypeError:  # if not a number that it is accepted. invalid point, anyways.
            self._coords = coords

    def _empty(self):
        """ revoke coordinates. This renders the Point invalid. """
        self._coords = EMPTY

    def _isempty(self):
        """ check if Point is empty (that is, invalid) """
        if self._coords == EMPTY:
            return True
        else:
            return False

    def distance(self, other):
        return distance(self, other)

    @cached_property
    def isvalid(self):
        """ checks for validity """
        if self._coords == EMPTY:
            # print('coords empty')
            return False
        if not isinstance(self.x, (int, float)):
            # print('coord x is not numeric')
            return False
        if not isinstance(self.y, (int, float)):
            # print('coord x is not numeric')
            return False
        if len(self.xy) != 2:
            # print('two coords are needed')
            return False
        return True

    def bounding_box(self):
        return (self.x, self.x), (self.y, self.y)

    def plot(self, show=False, style=None, annotate=False, txt=''):

        if _plotting_available:
            if style is None:
                self.plot_style = 'r.'
            else:
                self.plot_style = style
            plt.plot(self.x, self.y, self.plot_style, markersize=10)
            if annotate:
                ax = plt.gca()
                ax.annotate('%.2f, %.2f' % self.xy, xy=(self.x, self.y), xycoords='data', xytext=(1, 1),
                            textcoords='offset points', horizontalalignment='left', verticalalignment='bottom')

            if txt != '':
                ax = plt.gca()
                ax.annotate(txt, xy=(self.x, self.y), xycoords='data', xytext=(1, 1),
                            textcoords='offset points', horizontalalignment='left', verticalalignment='bottom')

            if show:
                plt.axis('tight')
                plt.axis('equal')
                plt.show()


class IntexPoint(Point):
    """
    This is a Point subclass that knows whick segment it is assigned to.
    If no segment is defined then it is practically a normal Point.
    """
    def __init__(self, args, segment=None):
        super(IntexPoint, self).__init__(args)
        self.segment = segment


def run():
    pass


if __name__ == '__main__':
    run()
