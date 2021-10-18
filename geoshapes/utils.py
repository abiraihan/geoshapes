# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 09:07:58 2021

@author: ABIR RAIHAN
"""

import numpy
import geopandas
import pandas
import shapely
from scipy.spatial import Voronoi


@classmethod
class Voronoi(
        cls,
        pointGeoms
        ):
    
    cls._pointGeoms = pointGeoms    
    
    if isinstance(cls._pointGeoms, str):
        pointData = geopandas.read_file(point_shape)
        if pointData.crs == None: pointData.crs = 'EPSG:4326'
    elif isinstance(cls._pointGeoms, geopandas.GeoDataFrame):
        pointData = cls._pointGeoms
        if pointData.crs == None: pointData.crs = 'EPSG:4326'
    else:
        raise TypeError("Expecting a string like object of path or geopandas.Geodataframe like object")
    
    convexHull = geopandas.GeoDataFrame(
        crs = 'EPSG:4326',
        geometry = [shapely.geometry.MultiPoint(
            [i for i in pointData.geometry]
            ).convex_hull]
        )
    
    listarray = []
    for pp in pointData.geometry:
        listarray.append([pp.x, pp.y])
    nparray = np.array(listarray)
    vor = Voronoi(nparray)
    lines = [
        shapely.geometry.LineString(vor.vertices[line])
        for line in vor.ridge_vertices
        if -1 not in line
        ]
    
    voronoiPoly = geopandas.GeoDataFrame()
    for poly in shapely.ops.polygonize(lines):
        singPoly = geopandas.GeoDataFrame(geometry=geopandas.GeoSeries(poly), crs = 'EPSG:4326')
        voronoiPoly = voronoiPoly.append(singPoly)
    
    voronoiPoly.reset_index(drop = True, inplace = True)
    return voronoiPoly, convexHull
    
    