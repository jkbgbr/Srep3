# -*- coding: utf-8 -*-

# try to import matplotlib, to enable visualisation
try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as PG
    from matplotlib.collections import PatchCollection as PC
    _plotting_available = True
    # plt.xkcd()
except ImportError:
    _plotting_available = False
    plt = None
    PG = None
    PC = None
    print('No plotting functions available, matplotlib is missing')

# As a global setting, if PolygonLines are NOT to be purged (see SR_polygon.purge()) set this to False.
# If PURGE == True, the polyline can be purged by explicitly calling purge() and is purged before transformations
# (rotate, move).
# Note, that tests are designed to pass with PURGE == True, PURGE == False is untested.
# If set to False, PolygonLines that have not been transformed MAY yield correct results. Transformed PolygonLines will
# most likely fail at some point during runtime.
# If set to False and care is taken to avoid defining coincident Points, the results will LIKELY be correct.
PURGE = True

# Help, Background material
# http://paulbourke.net/geometry/polygonmesh/

# triangle
# http://www.cs.cmu.edu/~quake/triangle.html
# http://dzhelil.info/triangle/index.html
# https://pypi.python.org/pypi/triangle/20170106
# http://block.arch.ethz.ch/blog/2016/06/using-jonathan-shewchuks-triangle-library-with-python/  INSTALL PROBLEM
