# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 15:56:00 2021

@author: ABIR RAIHAN
"""

import numpy as np
import geopandas
import pandas
import rasterio
from sklearn.neighbors import KNeighborsRegressor
from pykrige.ok import OrdinaryKriging
from sklearn.preprocessing import robust_scale
from sklearn.cluster import KMeans
from pysal.lib import weights
from pysal.explore import esda
from pysal.model import spreg

class interpolation:
    
    def __init__(self,gdf, attribute):
        
        
        self.gdf=gdf
        self.attribute = attribute
        self.x = gdf.geometry.x.values
        self.y = gdf.geometry.y.values
        self.z = gdf[attribute].values
        self.crs = gdf.crs
        self.xmax = gdf.geometry.x.max()
        self.xmin = gdf.geometry.x.min()
        self.ymax = gdf.geometry.y.max()
        self.ymin = gdf.geometry.y.min()
        self.extent = (self.xmin, self.xmax, self.ymin, self.ymax)
        self.res = (self.xmax - self.xmin) / 250
        self.ncol = int(np.ceil((self.xmax - self.xmin) / self.res)) # delx
        self.nrow = int(np.ceil((self.ymax - self.ymin) / self.res))# dely
        
    def points_to_grid(self):
        
        hrange = ((self.ymin,self.ymax),(self.xmin,self.xmax))
        zi, yi, xi = np.histogram2d(self.y, self.x, bins=(int(self.nrow), int(self.ncol)), weights=self.z, normed=False,range=hrange)
        counts, _, _ = np.histogram2d(self.y, self.x, bins=(int(self.nrow), int(self.ncol)),range=hrange)
        np.seterr(divide='ignore',invalid='ignore')
        zi = np.divide(zi,counts)
        np.seterr(divide=None,invalid=None)
        zi = np.ma.masked_invalid(zi)
        array = np.flipud(np.array(zi))
        return array
    
    def knn_2D(self, k=4, weights='uniform', algorithm='brute', p=2, maxrows = 9000):
        if len(self.gdf) > maxrows:
            raise ValueError('GeoDataFrame should not be larger than 9000 rows,\
                             knn is a slow algorithim and can be too much for your computer,\
                                 Change maxrows at own risk')  # shorthand for 'raise ValueError()'
        array = self.points_to_grid()
        X = []
        frow, fcol = np.where(np.isfinite(array))
        for i in range(len(frow)):
            X.append([frow[i], fcol[i]])
        y = array[frow, fcol]
        train_X, train_y = X, y
        knn = KNeighborsRegressor(n_neighbors=k, weights=weights, algorithm=algorithm, p=2)
        knn.fit(train_X, train_y)
        X_pred = []
        for r in range(int(self.nrow)):
            for c in range(int(self.ncol)):
                X_pred.append([r, c])
        y_pred = knn.predict(X_pred)
        karray = np.zeros((self.nrow, self.ncol))
        i = 0
        for r in range(int(self.nrow)):
            for c in range(int(self.ncol)):
                karray[r, c] = y_pred[i]
                i += 1
        return karray
    
    def OrdinaryKriging_2D(
            cls, n_closest_points=3,variogram_model='linear',
            verbose=False, coordinates_type='euclidean',
            backend='vectorized'):
        
        # if not pykrige_install:
        #     raise ValueError('Pykrige is not installed, try pip install pykrige')

        OK = OrdinaryKriging(cls.x,cls.y,cls.z, variogram_model=variogram_model, verbose=verbose,
                     enable_plotting=False, coordinates_type=coordinates_type)

        xpts = np.arange(cls.xmin + cls.res/2,cls.xmax+cls.res/2, cls.res)
        ypts = np.arange(cls.ymin + cls.res/2,cls.ymax+cls.res/2, cls.res)
        ypts = ypts[::-1]


        xp, yp = [],[]
        for yi in ypts:
            for xi in xpts:
                xp.append(xi)
                yp.append(yi)

        if n_closest_points is not None:
            backend = 'loop'
        krige_array, ss = OK.execute('points', xp, yp,n_closest_points=n_closest_points,backend=backend)

        krige_array = np.reshape(krige_array,(cls.nrow,cls.ncol))

        return krige_array
    
    def write_raster(self,array,path):
        if '.' not in path[-4:]:
            path += '.tiff'
        transform = rasterio.transform.from_origin(
            self.xmin, self.ymax,
            self.res, self.res
            )
        new_dataset = rasterio.open(
            path, 'w', driver='GTiff',
            height=array.shape[0], width=array.shape[1],
            count=1, dtype=array.dtype,
            crs=self.gdf.crs, transform=transform,
            nodata=np.nan)
        new_dataset.write(array, 1)
        new_dataset.close()
        return new_dataset
