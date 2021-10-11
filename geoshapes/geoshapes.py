# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 23:36:28 2021

@author: ABIR RAIHAN
"""
import os
import math
import numpy
import pandas
import shapely
import geopandas
import itertools
import pathlib

class splitShape:
    
    """
    To split Polygon geometry into different shape that required to create
    experimental plot / Trial design. It creates/re-creates experimental
    plot design to summerize data to geospatially enabled space for next-step
    mathematical modelling.
    """
    
    @classmethod
    def splitLatin(
            cls,
            geoms:shapely.geometry.Point,
            bufferLength:int
            ):
        
        """
        To split Polygon geometry into latin design to create a experimental plot.
        For reference : https://online.stat.psu.edu/stat503/lesson/4/4.3

        Parameters
        ----------
        geoms : shapely.geometry.Point
            DESCRIPTION. Single point shapely geometry
        bufferLength : int
            DESCRIPTION. Length of buffer or, side length

        Returns
        -------
        TYPE
            DESCRIPTION. List of shapely polygon geometry

        """
        
        cls._geoms = geoms
        cls._bufferLength = bufferLength
        
        xcords, ycords = cls._geoms.x, cls._geoms.y
        
        point_data = shapely.geometry.Point((xcords, ycords)).buffer(cls._bufferLength)
            
        s = point_data.simplify(0.5, preserve_topology=False)
        
        num_tiles = 3
        minx, miny, maxx, maxy = s.bounds
        dx = (maxx - minx) / num_tiles
        dy = (maxy - miny) / num_tiles
        
        lines = []
        for x in range(num_tiles + 1):
            lines.append(shapely.geometry.LineString([(minx + x * dx, miny), (minx + x * dx, maxy)]))
        for y in range(num_tiles + 1):
            lines.append(shapely.geometry.LineString([(minx, miny + y * dy), (maxx, miny + y * dy)]))
            
        return list(shapely.ops.polygonize(shapely.ops.unary_union(shapely.geometry.MultiLineString(lines))))
    
    @classmethod
    def splitCircle(
            cls,
            geoms:shapely.geometry.Point,
            circleRadius:float,
            incrementDegree:int,
            **kwargs:dict
            ):
        
        '''
        Parameters
        ----------
        To split a Polygon geometry into a Circle grid, which can be used in experimental design
        after selecting spatial location of interest.
        
        geoms : shapely.geometry.Point
            DESCRIPTION. single shapely Point geomtery
        circleRadius : float
            DESCRIPTION. circle radius in feet - float
        incrementDegree : int
            DESCRIPTION. degree increament step-wise
        **kwargs : dict
            DESCRIPTION. optional parametrs
        clipInterior : boolean
            Default is False. if True, returns intersected geomerty clipped out the innerWidth distance
        innerWidth : int
            Default is 1.
            Assign the number of feet that it should be intersected from the whole geometry from the center.
        getGeom : str
            'Inner', 'Outer' and 'Both': returns the geometry as assigned
            
        Returns
        -------
        list
            Returns a collection of shapely polygon geometry
        
        Purpose
        -------
            To get splited polgon for circular area. i.e. Experimental Plot / Trial Design
        
        Uses
        ----
            import string
            c = splitShape.splitCircle(point, 3188, 5, clipInterior = True, innerWidth = 2788, getGeom = 'Outer')
            ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:3857')
            ft['ids'] = range(len(ft))
            ft['Group']= ft.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
            ft.plot(cmap = 'tab20')

        '''
        
        cls._geoms = geoms
        cls._radius = circleRadius
        cls._degree = incrementDegree
        if cls._degree > 170:
            raise ValueError(f"Rotation degree not acceptable more than 170")
        
        data_args = dict(clipInterior = False,
                         innerWidth = 100,
                         getGeom = 'Both')
        
        for key, value in data_args.items():
            if key in kwargs:
                data_args[key] = kwargs[key]
        
        bufferLength = lambda x : x*0.3048
        
        rotateLine = lambda line, increment : [shapely.affinity.rotate(
            line,
            i,
            origin = line.coords[1])
            for i in numpy.arange(
                    0,
                    360,
                    float(increment)
                    )]
        
        geomClip = lambda singlePoly, bufferCut : [i.difference(
            bufferCut) for i in singlePoly]
        geomsBuffer = cls._geoms.buffer(bufferLength(cls._radius))
        
        UpperLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(
                geomsBuffer.bounds[0],
                geomsBuffer.bounds[3]),
             shapely.geometry.Point(
                 geomsBuffer.bounds[2],
                 geomsBuffer.bounds[3]
                 +(geomsBuffer.bounds[3]*bufferLength(
                                        cls._radius
                                        )))]).centroid
        
        incrementLines = rotateLine(
            shapely.geometry.LineString(
            [UpperLineCenter, cls._geoms]),
            cls._degree)
        
        cutBounds = []
        
        for i in range(len(incrementLines)):
            
            if i == len(incrementLines) - 1:
                Line1 = incrementLines[i].coords[0]
                Line2 = incrementLines[0].coords[0]
                auLine = incrementLines[i].coords[1]
                f = shapely.geometry.Polygon(
                    shapely.geometry.LinearRing(
                        [Line1, Line2, auLine]
                        ))
                cutBounds.append(f.intersection(geomsBuffer))
                
            else:
                Line1 = incrementLines[i].coords[0]
                Line2 = incrementLines[i+1].coords[0]
                auLine = incrementLines[i].coords[1]
                f = shapely.geometry.Polygon(
                    shapely.geometry.LinearRing(
                        [Line1, Line2, auLine]
                        ))
                cutBounds.append(f.intersection(geomsBuffer))
        
        if data_args['clipInterior'] == True:
            if data_args['innerWidth'] >= cls._radius:
                raise ValueError(f" 'InnerWidth' should be less than assigned radius for circle")
            if not data_args['getGeom'] in ['Inner', 'Outer', 'Both']:
                raise ValueError(f"{data_args['getGeom']} is not listed.\n\
                                 Only acceptable 'getGeom' parameter value is 'Inner', 'Outer' and 'Both'.")
            halfLength = (cls._radius + data_args['innerWidth'])/2
            clippedArea = geomClip(
                cutBounds,
                cls._geoms.buffer(
                    bufferLength(
                        data_args['innerWidth'])))
            bothArea = [list(shapely.ops.split(
                i, (cls._geoms.buffer(
                    bufferLength(halfLength)).boundary)))
                        for i in clippedArea]
            if data_args['getGeom'] == 'Both':
                return [k for i in bothArea for k in i]
            elif data_args['getGeom'] == 'Inner':
                return [i[1] for i in bothArea]
            elif data_args['getGeom'] == 'Outer':
                return [i[0] for i in bothArea]
        else:
            return cutBounds
    
    @classmethod
    def splitCircleSquare(
            cls,
            geoms:shapely.geometry.Point,
            circleRadius:float,
            rotation:int = 45
            ):
        
        """
        Parameters
        ----------
        geoms : shapely.geometry.Point
            DESCRIPTION. A single shapely Point geometry
        circleRadius : float
            DESCRIPTION. Insert the circle Radius in feet.
        rotation : int, optional
            DESCRIPTION. The default is 45.

        Returns
        -------
        TYPE
            DESCRIPTION. List of shapely polygon geometry

        """
        
        cls._geoms = geoms
        cls._radius = circleRadius
        cls._rotation = rotation
        
        bufferLength = lambda x : x*0.3048
        
        geomsBuffer = cls._geoms.buffer(
            bufferLength(cls._radius))
        
        xmin, ymin, xmax, ymax = geomsBuffer.bounds
        UpperLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmin,ymax),
             shapely.geometry.Point(xmax,ymax)]).centroid
        
        RightLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmax,ymax),
             shapely.geometry.Point(xmax,ymin)]).centroid
        
        LowerLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmax,ymin),
             shapely.geometry.Point(xmin,ymin)]).centroid
        
        LeftLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmin,ymin),
             shapely.geometry.Point(xmin,ymax)]).centroid
        
        polys = shapely.geometry.Polygon(
            shapely.geometry.LinearRing(
                [UpperLineCenter, RightLineCenter,
                 LowerLineCenter, LeftLineCenter]))
        
        polyRotation = shapely.affinity.rotate(
            polys,
            cls._rotation,
            origin = cls._geoms
            ).buffer(cls._radius/10000)
        
        polyLinestring = geomsBuffer.boundary.union(
            geomsBuffer.intersection(polyRotation).boundary)
        
        return list(shapely.ops.polygonize(polyLinestring))
    
    @classmethod
    def splitSquare(
            cls,
            geoms:shapely.geometry.Point,
            sideLength:float,
            rotation:int = 45,
            includeInterior:bool = True
            ):
        
        """
        Parameters
        ----------
        geoms : shapely.geometry.Point
            DESCRIPTION. A single shapely Point geometry 
        sideLength : float
            DESCRIPTION. Insert the sideLength of the square geometry in feet.
        rotation : int, optional
            DESCRIPTION. The default is 45.
        includeInterior : boolean, optional
            DESCRIPTION. The default is True. if 'False', returns polygon without the interior polygon shape
    
        Returns
        -------
        TYPE
            DESCRIPTION. List of shapely polygon geometry
    
        """
        
        cls._geoms = geoms
        cls._sideLength = sideLength
        cls._rotation = rotation
        cls._include = includeInterior
        
        bufferLength = lambda x : x*0.3048
        
        geomsBuffer = cls._geoms.buffer(
            bufferLength(
                cls._sideLength/numpy.sqrt(2)))
        xmin, ymin, xmax, ymax = geomsBuffer.bounds
        
        UpperLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmin,ymax),
             shapely.geometry.Point(xmax,ymax)]).centroid
        
        RightLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmax,ymax),
             shapely.geometry.Point(xmax,ymin)]).centroid
        
        LowerLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmax,ymin),
             shapely.geometry.Point(xmin,ymin)]).centroid
        
        LeftLineCenter = shapely.geometry.LineString(
            [shapely.geometry.Point(xmin,ymin),
             shapely.geometry.Point(xmin,ymax)]).centroid
        
        polys = shapely.affinity.rotate(
            shapely.geometry.Polygon(
                shapely.geometry.LinearRing(
                    [UpperLineCenter, RightLineCenter,
                     LowerLineCenter, LeftLineCenter])),
            45, origin = cls._geoms)
        
        polyRotation = shapely.affinity.rotate(
            polys, cls._rotation, origin = cls._geoms)
        
        geomsInterior = polyRotation.buffer(
            -bufferLength(
                cls._sideLength/3.618034))
        
        lines = shapely.geometry.LineString(
            [UpperLineCenter, LowerLineCenter])
        
        line1 = shapely.affinity.rotate(
            lines, cls._rotation, origin = cls._geoms)
        line2 = shapely.affinity.rotate(
            lines, cls._rotation + 90, origin = cls._geoms)
        
        split2 = [items
                  for subList in [list(i)
                                  for i in [shapely.ops.split(
                                          i,line2)
                                          for i in list(
                                                  shapely.ops.split(polyRotation,
                                                                    line1))]]
                  for items in subList]
        
        leftOver = [i.difference(geomsInterior)
                    for i in split2] + [geomsInterior]
        
        if cls._include == True:
            polyShapes = list(shapely.ops.polygonize(leftOver))
        else:
            polyShapes = [i.difference(geomsInterior)
                          for i in split2]
        
        return polyShapes
    
    @classmethod
    def splitGeom(
            cls,
            geoms:shapely.geometry.Polygon,
            splits:int,
            **kwargs:dict
            ):
        
        """
        Parameters
        ----------
        geoms : shapely.geometry.Polygon
            DESCRIPTION. A single shapely polygon geometry
        splits : int
            DESCRIPTION. number of splits required
        **kwargs : dict
            DESCRIPTION. 
        rotation : int
            DESCRIPTION. Rotation angle in degree, it will search the major axis for the polygon from 0 to assigned degree.
            Default is 30

        Returns
        -------
        TYPE
            DESCRIPTION. List of shapely polygon or multipolygon geometry
        
        Purpose
        -------
            Split a polygon by spliting number
            
        Usage
        ----
            sdf = geopandas.read_file("./geoDataFrame.shp")
            fl = shapely.geometry.box(*sdf.geometry[0].bounds).intersection(sdf.geometry[0])
            c = splitShape.splitGeom(geoms = fl, splits = 4, rotation = 45)
            ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:4326')
            ft.plot(cmap = 'Spectral')

        """
        
        cls._geoms = geoms
        cls._splits = splits
        
        options = dict(rotation= 30, singleGeom = False)
        
        for key, value in options.items():
            if key in kwargs:
                options[key] = kwargs[key]
        
        def getDegree(deg:int):
            if deg < 0:
                return 180 - (deg - ((deg//180)*180))
            elif 0<=deg<360:
                return deg - ((deg//180)*180)
            elif deg>=360:
                return deg - ((deg//180)*180)
            else:
                return 0
        
        if options['rotation'] != 30:
            
            options['rotation'] = getDegree(int(options['rotation']))
            print(f"  -- Found rotation angle {options['rotation']} degree")
        
        poly = shapely.geometry.Polygon(
            [shapely.geometry.Point((cls._geoms.bounds[0], cls._geoms.bounds[3])),
             shapely.geometry.Point((cls._geoms.bounds[0], cls._geoms.bounds[1])),
             shapely.geometry.Point((cls._geoms.bounds[2], cls._geoms.bounds[1])),
             shapely.geometry.Point((cls._geoms.bounds[2], cls._geoms.bounds[3])),
             shapely.geometry.Point((cls._geoms.bounds[0], cls._geoms.bounds[3]))]
            )
        print(f"  -- Analyzing suitable radians for the major axis at {options['rotation']} degree")
        bounddo =[shapely.affinity.rotate(poly, i)
                  for i in range(0, int(options['rotation']), 1)
                  if (shapely.affinity.rotate(poly, i).difference(cls._geoms)).area ==
                  min([(shapely.affinity.rotate(poly, i).difference(cls._geoms)).area
                       for i in range(0, int(options['rotation']), 1)])]
        print('  -- Analyzing Polygon Rotation is completed')

        x, y = list(bounddo[0].exterior.xy[0]), list(bounddo[0].exterior.xy[-1])
        coords = list(map(lambda x, y : (x, y), x, y))
        linefo = [shapely.geometry.LineString((shapely.geometry.Point(coords[i-1]), shapely.geometry.Point(coords[i])))
                  for i in range(len(coords)) if i!=0]
        print('  -- Major Axis Identified Successfully')
        try:
            rlist = [[i for i in linefo if i.length ==
                      max([i.length for i in linefo])][0].interpolate(i)
                     for i in [(max([i.length for i in linefo])/cls._splits)*i
                               for i in range(0, (cls._splits +1), 1)]]
            llist = [[i for i in linefo if i.length ==
                      max([i.length for i in linefo])][1].interpolate(i)
                     for i in [(max([i.length for i in linefo])/cls._splits)*i
                               for i in range(0, (cls._splits+1), 1)]][::-1]
        except IndexError:
            ld = geopandas.GeoDataFrame(
                geometry = [shapely.geometry.LineString((shapely.geometry.Point(coords[i-1]), shapely.geometry.Point(coords[i])))
                            for i in range(len(coords)) if i!=0])
            linefo1 = [ld.iloc[0].geometry, ld.iloc[1].geometry]
            rlist = [linefo1[0].interpolate(i)
                     for i in [(max([i.length
                                     for i in linefo1])/cls._splits)*i
                               for i in range(0, (cls._splits +1), 1)]]
            llist = [linefo1[1].interpolate(i)
                     for i in [(max([i.length
                                     for i in linefo1])/cls._splits)*i
                               for i in range(0, (cls._splits +1), 1)]][::-1]
        print('  -- Point location for polygon geometry identified successfully')
        splited_poly = geopandas.clip(geopandas.GeoDataFrame(
            crs = 'EPSG:4326',
            geometry = [shapely.geometry.Polygon([rlist[i], rlist[i+1], llist[i+1], llist[i]])
                        for i in range(len(rlist)-1)]),
            geopandas.GeoDataFrame(crs = 'EPSG:4326',
                             geometry = [cls._geoms]))
        print('  -- Geometry splited successfully and polygon defined accordingly')
        print(f"  -- Number of Splited Polygon before merge : {len(splited_poly)}")
        splited_poly.reset_index(drop = True, inplace = True)
        if options['singleGeom'] == True:
            splitedPolygon = [items
                              for sublist in [list(i)
                                              for i in splited_poly.geometry
                                              if i.geom_type == 'MultiPolygon']
                              for items in sublist]+ [i for i in splited_poly.geometry
                                                      if i.geom_type == 'Polygon']
            return splitedPolygon
        else:
            return splited_poly

class checkShape:
    
    @classmethod
    def checkPolygon(
            cls,
            data
            ):
        
        """
        Parameters
        ----------
        data : TYPE
            DESCRIPTION. geopandas GeoDataFrame or shapefile path

        Returns
        -------
        valiData : TYPE - geopandas geodataframe
            DESCRIPTION. return valid polygon geometry geodataframe
        
        Purpose
        -------
            DESCRIPTION. Chcek polygon geometry validity and returns only valid polygon geometry
            
        """
        cls.data = data
        if cls.data is None:
            raise ValueError(f"No data to read, data is {cls.data}")
            
        if isinstance(data, str):
            gdf = geopandas.read_file(cls.data)
            
        elif isinstance(data, geopandas.GeoDataFrame):
            gdf = cls.data
        else:
            pass
        
        kdf = gdf.set_geometry('geometry')
        df = kdf.geometry.explode()
        df.reset_index(drop = True, inplace = True)
        geoms = []
        for i in df.geometry:
            if not i.is_valid:
                make_valid = i.buffer(0)
                assert make_valid.geom_type == 'Polygon', f"Found Geometry Type {make_valid.geom_type} instead of Polygon"
                assert make_valid.is_valid
                geoms.append(make_valid)
                print(f' -- Invalid {make_valid.geom_type} geometry is fixed, after performing the fix geometry algorithm.')
            elif i.is_valid == True:
                geoms.append(i)
            else:
                print(f' -- Cannot Take Decision, avoided geometry')
                break
        valiData = geopandas.GeoDataFrame(
            index = [i for i in range(len(geoms))],
            crs = 'EPSG:4326',
            geometry = geoms)
        print(f"\n  -- Found No Issue : Performed Topology Check for the geometry")
        return valiData
    
    @classmethod
    def checkOverlaps(
            cls,
            geomPaths:str,
            pathType:str = 'directoryPath',
            tolerance:float = 0.005
            ):
        
        """
        Parameters
        ----------
        geomPaths : str
            DESCRIPTION. directory of file path location
        pathType : str, optional
            DESCRIPTION. The default is 'directoryPath'. Directory Location of shapefile.
            if 'filePath', single shapefile path location
        tolerance : float, optional
            DESCRIPTION. The default is 0.005 which is square feet.
            Insert the least square feet area that should returns in terms of overlap geometry.

        Raises
        ------
        ValueError
            DESCRIPTION. pathType should be either 'directoryPath' or 'filePath'.
        TypeError
            DESCRIPTION. Only accepted type of geometry is Polygon or MultiPolygon.
            If any other Type of geometry entered then returns TypeError

        Returns
        -------
        None. if not found any overlapped geometry from the data layer
        
        If found overlapped geometry into the data layer.
        
        returns '_checkOverlapGeometry' directory within the same directory/file path location with overlapped geometry.
        
        Purpose
        -------
            To identify and check the topology of polygon data layer.
        
        Uses
        ----
            If path is a file:
                checkShape.CheckOverlaps(filepath, 'filePath', 0.0005)
                
            If path is a  directory:
                checkShape.CheckOverlaps(directorypath)

        """
        
        cls._geomPaths = geomPaths
        cls._pathType = pathType
        cls._tolerance = tolerance
        
        if not cls._pathType in ['directoryPath', 'filePath']:
            raise ValueError("Only acceptable path type is either 'directoryPath' or 'filePath' keyword")
            
        if cls._pathType == 'directoryPath':
            geomfiles = [i.replace('/', '\\')
                         for i in [os.path.join(root, file)
                                   for root, dirs, files in os.walk(f"{cls._geomPaths}")
                                   for file in files if file.endswith(".shp")]]
            checkData = lambda file : pandas.concat(
                [geopandas.read_file(shp)
                 for shp in file]).pipe(geopandas.GeoDataFrame) if len(file) > 1 else geopandas.read_file(file[0])
            if len(geomfiles) != 0:
                geodata = checkData(geomfiles)
                geomData = geodata.to_crs('EPSG:4326')
            else:
                raise ValueError("Data location doesn't contain any shapefile, please check the 'pathType'/'geomPaths' that has been assigned ")
            
            pathGeoms = str(cls._geomPaths)
                
        elif cls._pathType == 'filePath':
            geodata = geopandas.read_file(cls._geomPaths)
            geomData = geodata.to_crs('EPSG:4326')
            pathGeoms = os.path.dirname(cls._geomPaths)
        else:
            raise ValueError("Data should be either directory or file path and explicitly attributed into the parameters")
        
        geoms = [items
                 for sublist in [list(i)
                                 for i in geomData.geometry
                                 if i.geom_type == 'MultiPolygon']
                 for items in sublist] + [i for i in geomData.geometry
                                          if i.geom_type == 'Polygon']
        
        if len(geoms) == 0:
            raise TypeError('Only accepted geometry is either Polygon or MultiPolygon')
            
        inter = shapely.ops.cascaded_union(
            [pair[0].intersection(pair[1])
             for pair in itertools.combinations(geoms, 2)]
            )
        collectGeoms = lambda geoms , leastArea: [i for i in list(geoms) if i.area*1e10*10.7639 >= float(leastArea)] if inter.geom_type == 'GeometryCollection' else [i for i in geoms if i.area*1e10*10.7639 >= float(leastArea)]
        
        geos = collectGeoms(inter, cls._tolerance)
        if len(geos) == 0:
            print(f'No overlapped polygon geometry found into the data layer')
        else:
            files = geopandas.GeoDataFrame(
                index = [i for i in range(len(geos))],
                crs = 'EPSG:4326',
                geometry = geos)
            path = os.path.join(
                f"{pathlib.Path(pathGeoms)}",
                "_checkOverlapGeometry")
            try: os.mkdir(path)
            except FileExistsError: pass
            files.to_file(f"{path}/overlapGeometry.shp")
            print(f' -- Please check the validity (Topology) of geometry data layer.\n\
                   -- Found overlapped geometry into the data layer.')
    
    @classmethod
    def checkGeometry(
            cls,
            fileObjects:geopandas.GeoDataFrame,
            geomType:str
            ):
        
        """
        Parameters
        ----------
        fileObjects : geopandas.GeoDataFrame
            DESCRIPTION. a geopandas geodatframe 
        geomType : str
            DESCRIPTION. geom_type name i.e. 'Polygon', 'Point', 'LineString'

        Returns
        -------
        bool
            DESCRIPTION. True if returns assigned geometry type else False

        """
        
        cls._fileObjects = fileObjects
        cls._geomType = geomType
        
        geomList = [i for i in cls._fileObjects.geometry if i.geom_type == cls._geomType]
        
        if len(geomList) > 0:
            return True
        else:
            return False

class gridShape:
    
    @classmethod
    def squareGrid(
            cls,
            boundary,
            length:int,
            cut:bool = False
            ):
        
        """        
        Parameters.
        ----------
        
        boundary : Boundary shapefile or geopandas GeoDataFrame type
            DESCRIPTION.
            Boundary shapefile or geopandas GeoDataFrame datatype is only allowed
            
        length : int
            DESCRIPTION.
            Length of the side of a square grid in meters
            
        Raises
        ------
        TypeError
            DESCRIPTION.
            If datatypes other than a shapefile path or Geopandas GeoDataFrame
            
        Returns
        -------
        TYPE : geopandas.GeoDataFrame
            DESCRIPTION.
            Square grid of the boundary
        
        Uses:
            
            bundary_shapefile_path = "./boundary_shapefile.shp"
            df = geopandas.read_file(bundary_shapefile_path)
            grid_side_length = 10
            clipped_gdf = createGrid.squareGrid(df, grid_side_length)
            or,
            clipped_gdf = createGrid.squareGrid(bundary_shapefile_path, grid_side_length)

        """
        cls._bounds = boundary
        cls._size = length
        cls._cut = cut
        
        if isinstance(cls._bounds, str):
            
            gdf = geopandas.read_file(cls._bounds)
            bnds = gdf.to_crs('EPSG:3857').geometry.buffer(cls._size).bounds
            
        elif isinstance(cls._bounds, geopandas.GeoDataFrame):
            
            gdf = cls._bounds
            bnds = gdf.to_crs('EPSG:3857').geometry.buffer(cls._size).bounds
            
        else:
            raise TypeError("Data types other than a shapefile path or Geopandas GeoDataFrame")
        
        x = numpy.linspace(bnds.minx[0], bnds.maxx[0], int(abs((bnds.minx[0]-bnds.maxx[0])/cls._size)))
        y = numpy.linspace(bnds.miny[0], bnds.maxy[0], int(abs((bnds.miny[0]-bnds.maxy[0])/cls._size)))
        
        hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(x[:-1], x[1:]) for yi in y]
        vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(y[:-1], y[1:]) for xi in x]
        grids = list(shapely.ops.polygonize(shapely.geometry.MultiLineString(hlines + vlines)))
        fg = geopandas.GeoDataFrame(
            index = [i for i in range(len(grids))],
            crs = 'EPSG:3857',
            geometry = grids)
        dfGrid = fg.to_crs('EPSG:4326')
        if cls._cut == True:
            return geopandas.clip(dfGrid, gdf)
        else:
            return dfGrid
    
    @classmethod
    def hexagonGrid(
            cls,
            gdf:geopandas.GeoDataFrame,
            spacing:int,
            cut:bool = False
            ):
        
        cls.gdf = gdf
        cls._spacing = spacing
        cls._cut = cut
        
        """
        Return a hexagon grid, based on the shape of *gdf* and on a *spacing* value (in
        units of *gdf*). If cut=False, the grid will not be intersected with *gdf*
        (i.e it makes a grid on the bounding-box of *gdf*).
        
        Parameters
        ----------
        gdf: GeoDataFrame
            The collection of polygons to be covered by the grid.
        spacing: Integer
            The dimension (will be used as the width) of the cells to create, in units of *gdf*.
        cut: Boolean, default True
            Cut the grid to fit the shape of *gdf* (cell partially covering it will be truncated). If False,
            the returned grid will fit the bounding box of *gdf*.
        Returns
        -------
        grid: GeoDataFrame
            A collection of polygons.
        """
        hSpacing = cls._spacing*1e-5
        vSpacing = cls._spacing*1e-5
        hOverlay = 0.0
        vOverlay = 0.0
    
        # To preserve symmetry, hspacing is fixed relative to vspacing
        xVertexLo = 0.288675134594813 * vSpacing
        xVertexHi = 0.577350269189626 * vSpacing
        hSpacing = xVertexLo + xVertexHi
        hOverlay = hSpacing
        halfVSpacing = vSpacing / 2.0
        
        xmin, ymin = [i.min() for i in gdf.bounds.T.values[:2]]
        xmax, ymax = [i.max() for i in gdf.bounds.T.values[2:]]
        cols = int(math.ceil((xmax - xmin) / hOverlay))
        rows = int(math.ceil((ymax - ymin) / (vSpacing - vOverlay)))
        
        geoms = []
        for col in range(cols):
            x1 = xmin + (col * hOverlay)
            x2 = x1 + (xVertexHi - xVertexLo)
            x3 = xmin + (col * hOverlay) + hSpacing
            x4 = x3 + (xVertexHi - xVertexLo)
    
            for row in range(rows):
                if (col % 2) == 0:
                    y1 = ymax + (row * vOverlay) - (((row * 2) + 0) * halfVSpacing)
                    y2 = ymax + (row * vOverlay) - (((row * 2) + 1) * halfVSpacing)
                    y3 = ymax + (row * vOverlay) - (((row * 2) + 2) * halfVSpacing)
                else:
                    y1 = ymax + (row * vOverlay) - (((row * 2) + 1) * halfVSpacing)
                    y2 = ymax + (row * vOverlay) - (((row * 2) + 2) * halfVSpacing)
                    y3 = ymax + (row * vOverlay) - (((row * 2) + 3) * halfVSpacing)
    
                geoms.append((
                    (x1, y2),
                    (x2, y1), (x3, y1), (x4, y2), (x3, y3), (x2, y3),
                    (x1, y2)
                ))
                
        if cls._cut == True:
            hexData = geopandas.GeoDataFrame(
                index=[i for i in range(len(geoms))],
                geometry=pandas.Series(geoms).apply(lambda x: shapely.geometry.Polygon(x)),
                crs=gdf.crs)
            return geopandas.overlay(
                gdf,
                hexData,
                how='intersection')
        else:
            return geopandas.GeoDataFrame(
                index=[i for i in range(len(geoms))],
                geometry=pandas.Series(geoms).apply(lambda x: shapely.geometry.Polygon(x)),
                crs=gdf.crs)
    
    @classmethod
    def sanitizeGrid(
            cls,
            gdf,
            tolerance = 11
            ):
        """

        Parameters
        ----------
        gdf : TYPE
            DESCRIPTION. A collection of polygon geometry as a geodataframe object
        tolerance : TYPE, optional
            DESCRIPTION. The default is 11.

        Raises
        ------
        ValueError
            DESCRIPTION. Tolerance value can not be less than 11 (1.1 centimeter).

        Returns
        -------
        sanitized : geopandas.GeoDataFrame
            DESCRIPTION. A collection of polygon geometry as a geodataframe object

        """
        cls._gdf = gdf
        cls._tolerance = tolerance
        if cls._tolerance < 11:
            raise ValueError(f" -- Tolerance can not be less than 11 (1.1 centimeter). \
                             Please assign tolerance value more than 11")
        
        if cls._gdf.crs == 'EPSG:4326':
        
            geoms = [i.buffer(1e-6).buffer(-cls._tolerance/1e9)
                     for i in [items for sublist in [list(i)
                                                     for i in cls._gdf.geometry if i.geom_type == 'MultiPolygon']
                               for items in sublist]
                     + [i for i in cls._gdf.geometry if i.geom_type == 'Polygon']]
        
            sanitized = geopandas.GeoDataFrame(
                index = [i for i in range(len(geoms))],
                crs = 'EPSG:4326',
                geometry = geoms)
        else:
            toGdf = cls._gdf.to_crs('EPSG:4326')
            geoms = [i.buffer(1e-6).buffer(-cls._tolerance/1e9)
                     for i in [items for sublist in [list(i)
                                                     for i in toGdf.geometry if i.geom_type == 'MultiPolygon']
                               for items in sublist]
                     + [i for i in toGdf.geometry if i.geom_type == 'Polygon']]    
        
            sanitized = geopandas.GeoDataFrame(
                index = [i for i in range(len(geoms))],
                crs = 'EPSG:4326',
                geometry = geoms)
        return sanitized

class mergeShape:
    
    @classmethod
    def mergePolygon(
            cls,
            geomData,
            mergedPoly:int
            ):
        """

        Parameters
        ----------
        geomData : list/geopandas.GeoDataFrame
             - List of shapely polygon geometry or geopandas GeoDataFrame.
        mergedPoly : int
            Least Number of Merged geometry.

        Raises
        ------
        TypeError
            - If data is not a list of geometry / A geopandas GeoDataFrame

        Returns
        -------
        saveGeom : geopandas.GeoDataFrame
            A geopandas GeoDataFrame

        """
        
        cls._geomData = geomData
        cls._mergedPoly = mergedPoly
        
        if isinstance(cls._geomData, list):
            allPoly = geopandas.geoDataFrame(geometry = cls._geomData)
            allPolyAreas = shapely.ops.unary_union(cls._geomData).area
        elif isinstance(cls._geomData, geopandas.GeoDataFrame):
            allPoly = cls._geomData
            allPolyAreas = shapely.ops.unary_union([i for i in cls._geomData.geometry]).area
        else:
            raise TypeError("Data types other than a list of geometry or Geopandas GeoDataFrame")
        
        equalArea = allPolyAreas/int(cls._mergedPoly)
        geoms , geom_area = [], 0
        saveGeom = geopandas.GeoDataFrame()
        for index in range(len(allPoly)):
            if index == len(allPoly) - 1:
                geoms.append(allPoly['geometry'][index])
                df = geopandas.GeoDataFrame(geometry = geoms)
                bounds = geopandas.GeoDataFrame(
                    crs = 'EPSG:4326',
                    geometry=geopandas.GeoSeries(
                        shapely.ops.unary_union(df.geometry).intersection(
                        shapely.ops.unary_union(df.geometry))))
                saveGeom = saveGeom.append(bounds)
            elif equalArea>= geom_area:
                geom_area += allPoly['geometry'][index].area
                geoms.append(allPoly['geometry'][index])
            else:
                geoms.append(allPoly['geometry'][index])
                df = geopandas.GeoDataFrame(geometry = geoms)
                bounds = geopandas.GeoDataFrame(
                    crs = 'EPSG:4326',
                    geometry=geopandas.GeoSeries(
                        shapely.ops.unary_union(df.geometry).intersection(
                        shapely.ops.unary_union(df.geometry))))
                saveGeom = saveGeom.append(bounds)
                geoms, geom_area = [], 0
        saveGeom.reset_index(drop = True, inplace = True)
        return saveGeom
    
    @classmethod
    def mergeSilevers(
            cls,
            geomData,
            splitGeoms:int
            ):
        """

        Parameters
        ----------
        geomData : list
            - A collection of geometry.
        splitGeoms : int
            - Least number of geometry that geomData will be splited

        Raises
        ------
        TypeError
            - If data is not a list of geometry / A geopandas GeoDataFrame

        Returns
        -------
        mergedGeoms : list
            - A list of shapely polygon collection

        """
        
        cls._geomData = geomData
        cls._splitGeoms = splitGeoms
        
        if isinstance(cls._geomData, list):
            listPoly = cls._geomData
        elif isinstance(cls._geomData, geopandas.GeoDataFrame):
            listPoly = [i for i in cls._geomData.geometry]
        else:
            raise TypeError("Data is not a list of geometry / A geopandas GeoDataFrame")
        
        areas = min(sorted([i.area for i in listPoly])[-cls._splitGeoms:])
        
        msdDict = {i.area : i for i in listPoly}
        
        originalsDetected = {i : j for i, j in msdDict.items() if i>=areas}
        sileversDetected = {i : j for i, j in msdDict.items() if i<areas}
        
        if len(sileversDetected) != 1:
            sileversDetected = {i.area : i for i in [shapely.ops.unary_union([i[1] for i in sileversDetected.items()])]}
        
        mergedGeoms = []
        for _, j in sileversDetected.items():
            for _, l in originalsDetected.items():
                if j.touches(l) == True:
                    geoms = shapely.ops.unary_union([j, l])
                    mergedGeoms.append(geoms)
                else:
                    mergedGeoms.append(l)
        
        return mergedGeoms
    
    @classmethod
    def mergeOverlaps(
            cls,
            sourceGeoms:shapely.geometry.Polygon,
            processedGeoms:list
            ):
        """

        Parameters
        ----------
        sourceGeoms : shapely.geometry.Polygon
            - A single shapely Polygon geometry.
        processedGeoms : list
            - A collection of shapely Polygon geometry.

        Returns
        -------
        list
            - A collection of shapely Polygon geometry.

        """
        
        cls._sourceGeoms = sourceGeoms
        cls._processedGeoms = processedGeoms
        
        sourceArea = cls._sourceGeoms.area*1e10 + 1
        print(f"  --| Source Geometry Area is : {round(sourceArea, 2)} square meter")
        processArea = sum([i.area*1e10 for i in cls._processedGeoms])
        print(f"  --| Processed Geometry Area is : {round(processArea, 2)} square meter")
        if sourceArea < processArea:
            gty = [i for i in cls._processedGeoms]
            
            lf = []
            for i in itertools.combinations(gty, 2):
                if i[0].overlaps(i[1]):
                    geo = shapely.ops.unary_union([i[0], i[1]])
                    if i[0].area>i[1].area:
                        lf.append(geo.difference(i[1]))
                        lf.append(i[1])
                    else:
                        lf.append(geo.difference(i[0]))
                        lf.append(i[0])
            return [i for i in gty for j in lf if not i.intersects(j)] + lf
        else:
            return cls._processedGeoms

class selectGeom:
    
    @classmethod
    def getGeom(
            cls,
            geom_data
            ):
        """

        Parameters
        ----------
        geom_data : shapely.geometry.Polygon
            - A shapely Polygon geometry

        Returns
        -------
        geopandas.GeoDataFrame
            - A geopandas GeoDataFrame with ratio information of max and min side length.

        """
        
        cls.geom_data = geom_data
        
        poly = shapely.geometry.Polygon(
            [shapely.geometry.Point((cls.geom_data.bounds[0], cls.geom_data.bounds[3])),
             shapely.geometry.Point((cls.geom_data.bounds[0], cls.geom_data.bounds[1])),
             shapely.geometry.Point((cls.geom_data.bounds[2], cls.geom_data.bounds[1])),
             shapely.geometry.Point((cls.geom_data.bounds[2], cls.geom_data.bounds[3])),
             shapely.geometry.Point((cls.geom_data.bounds[0], cls.geom_data.bounds[3]))]
            )
        
        x, y = list(poly.exterior.xy[0]), list(poly.exterior.xy[-1])
        linefo = [shapely.geometry.LineString((shapely.geometry.Point(list(map(lambda x, y : (x, y), x, y))[i-1]), shapely.geometry.Point(list(map(lambda x, y : (x, y), x, y))[i])))
                  for i in range(len(list(map(lambda x, y : (x, y), x, y)))) if i!=0]
        mk = [i for i in linefo if i.length == max([linefo[i-1].length for i in range(len(linefo))])][0]
        lk = [i for i in linefo if i.length == min([linefo[i-1].length for i in range(len(linefo))])][0]
        df = geopandas.GeoDataFrame(crs = 'EPSG:4326', geometry = [cls.geom_data])
        df['max_length'],  df['min_length'] = mk.length*1e5, lk.length*1e5
        df['ratio'] = df.apply(lambda x : (x['max_length']/x['min_length']), axis = 1)
        return df
    
    @classmethod
    def getPolygon(
            cls,
            data:geopandas.GeoDataFrame,
            shape:str,
            sample:int
            ):
        """

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            - A geopandas GeoDataFrame for a collection of Polygon.
        shape : str
            - shape Name; Accepted keyword is either ['Square', 'Rectangle', 'Square-Rectangle']
        sample : int
            - Number of sample required

        Raises
        ------
        AttributeError
            - If Accepted keyword ['Square', 'Rectangle', 'Square-Rectangle'] is not mentioned.

        Returns
        -------
        geopandas.GeoDataFrame / shapeError
            - A geopandas GeoDataFrame

        """
        
        cls.data = data
        cls.shape = shape
        cls.sample = sample
        
        len_data = geopandas.GeoDataFrame()
        for index, row in cls.data.iterrows():
            ldata = cls.getGeom(row['geometry'])
            len_data = len_data.append(ldata)
        len_data.reset_index(drop = True, inplace = True)
        # if 'index' in len_data.columns.unique().tolist():del len_data['index']
        # len_data.crs = 'EPSG:4326'
        if cls.shape == 'Rectangle':
            return len_data.iloc[len_data['ratio'].argsort()[-cls.sample:]], len_data.iloc[~len_data['ratio'].argsort()[-cls.sample:]]
        elif cls.shape == 'Square':
            return len_data.iloc[len_data['ratio'].argsort()[:cls.sample]], len_data.iloc[~len_data['ratio'].argsort()[:cls.sample]]
        elif cls.shape == 'Square-Rectangle':
            return len_data.iloc[len_data['ratio'].argsort()[int((len_data['ratio'].max())/2):int((len_data['ratio'].max())/2)+cls.sample]], len_data.iloc[~len_data['ratio'].argsort()[int((len_data['ratio'].max())/2):int((len_data['ratio'].max())/2)+cls.sample]]
        else:
            raise AttributeError("Not accepted, only acceptable item are 'Rectangle, Square, Square-Rectangle'")
    
    @classmethod
    def getInteriorGeoms(
            cls,
            geomData,
            interiorBufferFeet:int = 1
            ):
        """
        Parameters
        ----------
        geomData : list/geopandas.GeoDataFrame
            - A geopandas GeoDataFrame or a list of shapely polygon geometry.
        interiorBufferFeet : int, optional
            - Buffer distance from polygon boundary. The default is 1.

        Raises
        ------
        TypeError
            - if data other than GeoDataFrame or list of shapely geometry.

        Returns
        -------
        list/None
            - Interior geometry of a list of shapely geometry or None
        list
            - Exterior geometry of a list of shapely geometry 

        """
        
        cls._geomData = geomData
        cls._interiorBufferFeet = interiorBufferFeet
        
        intBuffer = (1 + int(cls._interiorBufferFeet))*3.048e-6
        
        if isinstance(cls._geomData, geopandas.GeoDataFrame):
            data = [i for i in cls._geomData.geometry]
        elif isinstance(cls._geomData, list):
            data = cls._geomData
        else:
            raise TypeError("Expecting a list of geometry collection or geopandas GeoDataFrame.")
        
        bounds = shapely.ops.unary_union(data).buffer(3.048e-6).buffer(-float(intBuffer))
        
        exterior = [i for i in data if bounds.boundary.intersects(i)]
        interior = [i for i in data if not bounds.boundary.intersects(i)]
        
        testBounds = shapely.ops.unary_union(data).buffer(3.048e-6).buffer(-float(6.096e-06))
        
        validOutputInterior = [i for i in interior if testBounds.boundary.intersects(i)]
        
        if len(validOutputInterior) == 0 and len(interior) != 0:
            print("  -- Found both Interior & Exterior geometry")
            return interior, exterior
        elif len(validOutputInterior) == 0 and len(interior) == 0:
            print("  -- No Interior, Only Exterior geometry filtered")
            return None, exterior
        else:
            print("  -- Not determined, Source geometry returned as Exterior geometry after processing")
            return None, data

class utils:
    
    @classmethod
    def singlePolygon(
            cls,
            geomsData
            ):
        """

        Parameters
        ----------
        geomsData : list/geopandas.GeoDataFrame
            - A list of shapely polygon geometry / A geopandas GeoDataFrame

        Raises
        ------
        TypeError
            - If data other than GeoDataFrame or list of shapely geometry.

        Returns
        -------
        list
            - A list of shapely Polygon geometry

        """
        
        cls._geomsData = geomsData
        
        if isinstance(cls._geomsData, list):
            data = cls._geomsData
        elif isinstance(cls._geomsData, geopandas.GeoDataFrame):
            data = [i for i in cls._geomsData.geometry]
        else:
            raise TypeError("Expecting a list of geometry \ collection or geopandas GeoDataFrame.")
        return [items for sublist in [list(i)
                                      for i in data
                                      if i.geom_type == 'MultiPolygon']
                for items in sublist
                ] + [i for i in data if i.geom_type == 'Polygon']
    
    @classmethod
    def mergeData(
            cls,
            dataPath:str
            ):
        """

        Parameters
        ----------
        dataPath : str
            - data dircetory location for geospatially enabled data

        Raises
        ------
        AttributeError
            - If no geospatially enabled data found into the assigned data directory

        Returns
        -------
        geopandas.GeoDataFrame
            - A geopandas GeoDataFrame

        """
        
        cls._dataPath = dataPath
        geomfiles = [i.replace('/', '\\')
                     for i in [os.path.join(root, file)
                               for root, dirs, files in os.walk(f"{cls._dataPath}")
                               for file in files
                               if file.endswith(".shp")]]
        checkData = lambda file : pandas.concat([
            geopandas.read_file(shp)
            for shp in file]).pipe(geopandas.GeoDataFrame)if len(file)>1 else geopandas.read_file(file[0])
        if len(geomfiles) != 0:
            return checkData(geomfiles)
        else:
            raise AttributeError("Data location doesn't contain any shapefile")
    
    @classmethod
    def delDuplicates(
            cls,
            gdf:geopandas.GeoDataFrame
            ):
        """
        Removes duplicate geometry along with its values from the geopandas GeoDataFrame

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            - A geopandas GeoDataFrame

        Returns
        -------
        dupGdf:geopandas.GeoDataFrame

        """
        cls._gdf = gdf
        cls._gdf["geometry"] = cls._gdf["geometry"].apply(lambda geom: geom.wkb)
        dupGdf = cls._gdf.drop_duplicates(["geometry"])
        dupGdf["geometry"] = dupGdf["geometry"].apply(lambda geom: shapely.wkb.loads(geom))
        return dupGdf