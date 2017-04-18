# import triangle
# import triangle.plot as plot
# import matplotlib.pyplot as plt
#
# # delaunay
# d0 = triangle.get_data('dots')
# d1 = triangle.triangulate(d0)
# triangle.plot.compare(plt, d0, d1)
# plt.show()
#
# # voronoi diagram
# pts = triangle.get_data('dots')['vertices']
# ax1 = plt.subplot(121, aspect='equal')
# triangle.plot.plot(ax1, vertices=pts)
# lim = ax1.axis()
# points, edges, ray_origin, ray_direct = triangle.voronoi(pts)
# d = dict(vertices=points, edges=edges, ray_origins=ray_origin, ray_directions=ray_direct)
# ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
# triangle.plot.plot(ax2, **d)
# ax2.axis(lim)
# plt.show()
#
# # Planar Straight Line Graph (PSLG)
# A = triangle.get_data('A')
# ax = plt.axes()
# plot.plot(ax, **A)
# plt.show()
#
#
# # constrained delaunay
# ax = plt.axes()
# A = triangle.get_data('A')
# t = triangle.triangulate(A, 'p')
# plot.plot(ax, **t)
# plt.show()
#
# # constrained conforming delaunay
# A = triangle.get_data('A')
# t = triangle.triangulate(A, 'pq30')
# plot.plot(plt.axes(), **t)
# plt.show()
