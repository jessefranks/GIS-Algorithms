# Date: 4/18/2017
# Final Project
# Name: Jesse Franks

import matplotlib.pyplot as plt
import numpy as np
import random
import sys
sys.path.append('C:\\Users\\Jesse\\Downloads\\!classes\\GEOG 5222\\lib\\geom')
from point import *
from extent import *
from kdtree1 import *
from kdtree2b import *
from kfunction import *
from osgeo import ogr

# Get shapfiles for processing.
driver = ogr.GetDriverByName("ESRI Shapefile")

fname = 'C:\\Users\\Jesse\\Downloads\\!classes\\GEOG 5222\\final project\\airport_final_2.shp' # USA Airport locations(point).
airPorts = driver.Open(fname, 0)   # 0 for read only, 1 for writable
# Imported shapfile has been projected to USA_Contiguous_Equidistant_Conic.


airPortLayer = airPorts.GetLayer(0) # Get layer with airport features. 

# make airport features into point class.
airportPointClass = []
for i in range(airPortLayer.GetFeatureCount()):
     f = airPortLayer.GetFeature(i)
     geom = f.GetGeometryRef()
     airportPointClass.append(Point(geom.GetPoint(0)[0], geom.GetPoint(0)[1])) # Point class



# Draw airports
def draw_points(points, area):
    fig, ax = plt.subplots()
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    ax.scatter(xs, ys, edgecolor='none', facecolor='green', alpha=0.5)
    plt.xlabel('meters')
    plt.title('Airports in South-East Region of US')
    plt.grid()
    plt.xlim(extent[0], extent[1])
    plt.ylim(extent[2], extent[3])
    ax.set_aspect(1)
    plt.show()
    

# Set area
extent = airPortLayer.GetExtent()
area = Extent(extent[0], extent[1], extent[2], extent[3]) # get area of airports from geomotry.

# Plot airports
draw_points(airportPointClass, area)
# Number of airports
n = len(airportPointClass)

# k-fucntion
def get_kfunction_values(points, area):
    n = len(points)
    density = float(n)/area.area()
    t = kdtree2(points)
    d = min([area.xmax-area.xmin, area.ymax-area.ymin])*2.0/3/10.0 # a 10th of 2/3 of the length
    ds = [ d*(i+1) for i in range(10)]
    lds = [0 for d in ds]
    for i, d in enumerate(ds):
        for p in points:
            ld = kfunc(t, p, d, density)[1]
            lds[i] += ld
    return ds, [ld/n for ld in lds]


#Monte Carlo test
def kfunc_monte_carlo(n, area, radii, density, rounds=100):
    """
    Input
      n:            number of points, not used
      area:         Extent object defining the area
      radii:        list containing a set of radii of circles
      density:      density of point events in the area
      rounds:       number of simulations
    Return
      percentiles:  a list of 2.5th and 97.5th percentiles
                    for each d in radii
    """
    alllds = []
    for test in range(rounds):
        N = np.random.poisson(n)
        x = list(np.random.uniform(area.xmin, area.xmax, N))
        y = list(np.random.uniform(area.ymin, area.ymax, N))
        points = [Point(x[i], y[i]) for i in range(N)]
        t = kdtree2(points)
        lds = [0 for d in radii]
        for i, d in enumerate(radii):
            for p in points:
                ld = kfunc(t, p, d, density)[1]
                lds[i] += ld
        lds = [ld/N for ld in lds]
        alllds.append(lds)
    alllds = np.array(alllds)
    percentiles = []
    for i in range(len(radii)):
        percentiles.append([np.percentile(alllds[:,i], 2.5),
                            np.percentile(alllds[:,i], 97.5)])
    return percentiles

# Call functions
ds, lds1 = get_kfunction_values(airportPointClass,area)
density = float(n)/area.area()
percentiles = kfunc_monte_carlo(n, area, ds, density)


#plot k-function
fig, ax = plt.subplots()
plt.plot(ds, [p[1] for p in percentiles], color='grey', label='97.5% Envelope')
plt.plot(ds, lds1, color='blue', label='Airports')
plt.plot(ds, [p[0] for p in percentiles], color='grey', label='2.5% Envelope')
plt.plot([0, plt.xlim()[1]], [0, plt.xlim()[1]], color='red', label='L(d)=d')

plt.xlabel('Radius ($d$)')
plt.ylabel('$L(d)$')
plt.title('K Function ($\lambda$ = {0:5.2f})'.format(density))
plt.legend(loc='right', bbox_to_anchor=(1.55, 0.5))

plt.xlabel('Radius ($d$)')
plt.ylabel('$L(d)$')
ax.set_aspect(1)
plt.show()




