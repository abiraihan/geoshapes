# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 07:32:29 2021

@author: ABIR RAIHAN
"""

import geopandas
import pandas
import shapely
import fiona


class readData:
    
    @classmethod
    def readExcel(
            cls,
            excelFilePath:str,
            sheetNames:str,
            **kwargs
            ):
        
        cls._excelFilePath = excelFilePath
        cls._sheetNames = sheetNames
        
        data_args = dict(
            shapeGeometry = 'Point',
            latAttributes = None,
            longAttributes = None,
            polyAttributes = None
            )
        
        for key, value in data_args.items():
            if key in kwargs: data_args[key] = kwargs[key]
        
        if not data_args['shapeGeometry'] in ['Point', 'Polygon']:
            raise ValueError("Only acceptable geometry can be created is Point or, Polygon")
        
        if data_args['latAttributes'] is None and data_args['shapeGeometry'] == 'Point':
            raise AttributeError(f"Please specify the Latitude attribute for {data_args['shapeGeometry']} geometry")
        
        if data_args['longAttributes'] is None and data_args['shapeGeometry'] == 'Point':
            raise AttributeError(f"Please specify the Longitude attribute for {data_args['shapeGeometry']} geometry")
        
        if data_args['polyAttributes'] is None and data_args['shapeGeometry'] == 'Polygon':
            raise AttributeError(f"Please specify the Polygon geometry attribute for {data_args['shapeGeometry']} geometry")
        
        dataFrame = pandas.read_excel(
            cls._excelFilePath,
            sheet_name = cls._sheetNames
            )
        if data_args['shapeGeometry'] == 'Point':
            geoData = geopandas.GeoDataFrame(
                dataFrame,
                geometry = geopandas.points_from_xy(
                    dataFrame[data_args['latAttributes']],
                    dataFrame[data_args['longAttributes']]),
                crs="EPSG:4326"
                )
        elif data_args['shapeGeometry'] == 'Polygon':
            geoData = geopandas.GeoDataFrame(
                dataFrame,
                geometry = dataFrame[data_args['polyAttributes']],
                crs="EPSG:4326"
                )
        else:
            raise ValueError("Only acceptable geometry can be created is Point or, Polygon")
        return geoData