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
import exception

class splitShape:
    
    @classmethod
    def splitLatin(
            cls,
            geoms:shapely.geometry.Point,
            bufferLength:int
            ):
        
        """

        Parameters
        ----------
        cls : TYPE
            DESCRIPTION. class method
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
        cls : TYPE
            DESCRIPTION. class method
        geoms : shapely.geometry.Point
            DESCRIPTION. single shapely Point geomtery
        circleRadius : float
            DESCRIPTION. circle radius in feet - float
        incrementDegree : int
            DESCRIPTION. degree increament step-wise
        **kwargs : dict
            DESCRIPTION. optional parametrs
        clipInterior : boolean
            Default is False. if True, returns intersected geomerty
        innerWidth : int
            Default is 1.
            Assign the number of feet that it should be intersected from the whole geometry.
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
        cls : TYPE
            DESCRIPTION. class method
        geoms : shapely.geometry.Point
            DESCRIPTION. single shapely Point geometry
        circleRadius : float
            DESCRIPTION. Insert the circleRadius in feet.
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
        cls : TYPE
            DESCRIPTION.
        geoms : shapely.geometry.Point
            DESCRIPTION. single shapely Point geometry 
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
        cls : TYPE
            DESCRIPTION. class method
        geoms : shapely.geometry.Polygon
            DESCRIPTION. single shapely polygon geometry
        splits : int
            DESCRIPTION. number of splits required
        **kwargs : dict
            DESCRIPTION. 
        rotation : int
            DESCRIPTION. Rotation angle in degree, insert the degree that required.
            Default is 30

        Returns
        -------
        TYPE
            DESCRIPTION. List of shapely polygon or multipolygon geometry
        
        Purpose
        -------
            Split any polygon by spliting number
            
        Usage
        ----
            sdf = geopandas.read_file("./geoDataFrame.shp")
            fl = shapely.geometry.box(*sdf.geometry[0].bounds).intersection(sdf.geometry[0])
            c = splitShape.splitGeom(fl, 4)
            ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:4326')
            ft.plot(cmap = 'Spectral')

        """
        
        cls._geoms = geoms
        cls._splits = splits
        
        options = dict(rotation= 30)
        
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
        
        splitedPolygon = [items
                          for sublist in [list(i)
                                          for i in splited_poly.geometry
                                          if i.geom_type == 'MultiPolygon']
                          for items in sublist]+ [i for i in splited_poly.geometry
                                                  if i.geom_type == 'Polygon']
        
        return splitedPolygon

class checkShape:
    
    @classmethod
    def checkPolygon(
            cls,
            data
            ):
        
        """
        Parameters
        ----------
        cls : TYPE
            DESCRIPTION. Class method type
        data : TYPE
            DESCRIPTION. geopandas geodataframe or shapefile path

        Returns
        -------
        valiData : TYPE - geopandas geodataframe
            DESCRIPTION. return valid polygon geometry geodataframe
        
        Purpose
        -------
            DESCRIPTION. Chcek polygon geometry validity and returns only valid polygon geometry
            
        """
        cls.data = data
        if isinstance(data, str):
            gdf = geopandas.read_file(cls.data)
            
        elif isinstance(data, geopandas.GeoDataFrame):
            gdf = cls.data
        else:
            pass
        gdf = gdf.set_geometry('geometry')
        df = gdf.geometry.explode()
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
        print(f"\n  --  Found No Issue : Performed Topology Check for the geometry")
        return valiData
    
    @classmethod
    def checkOverlaps(
            cls,
            geomPaths,
            pathType:str = 'directoryPath',
            tolerance:float = 0.005
            ):
        
        """
        Parameters
        ----------
        cls : TYPE
            DESCRIPTION. Class Method Type
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
                raise ValueError("Data location doesn't contain any shapefile, please check the 'pathType' that has been assigned ")
            
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
        cls : TYPE
            DESCRIPTION. Class Method
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
        InputError
            DESCRIPTION.
            If datatypes other than a shapefile path or Geopandas GeoDataFrame
            
        Returns
        -------
        TYPE : geopandas.GeoDataFrame
            DESCRIPTION.
            Square grid of the boundary
        
        Uses:
            
            bundary_shapefile_path = "./boundary_shapefile.shp"
            
            df = gpd.read_file(bundary_shapefile_path)
            
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
            raise exception.InputError(type(cls._bounds), "Shape geometry")
        
        x = numpy.linspace(bnds.minx[0], bnds.maxx[0], int(abs((bnds.minx[0]-bnds.maxx[0])/cls._size)))
        y = numpy.linspace(bnds.miny[0], bnds.maxy[0], int(abs((bnds.miny[0]-bnds.maxy[0])/cls._size)))
        
        hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(x[:-1], x[1:]) for yi in y]
        vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(y[:-1], y[1:]) for xi in x]
        grids = list(shapely.ops.polygonize(shapely.geometry.MultiLineString(hlines + vlines)))
        fg = geopandas.GeoDataFrame(crs = 'EPSG:3857', geometry = grids)
        dfGrid = fg.to_crs('EPSG:4326')
        if cls._cut == True:
            return geopandas.clip(dfGrid, gdf)
        else:
            return dfGrid
    
    @classmethod
    def hexagonGrid(
            self,
            gdf:geopandas.GeoDataFrame,
            spacing:int,
            cut:bool = False
            ):
        
        self.gdf = gdf
        self._spacing = spacing
        self._cut = cut
        
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
        hSpacing = self._spacing*1e-5
        vSpacing = self._spacing*1e-5
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
            # (column + 1) and (row + 1) calculation is used to maintain
            # topology between adjacent shapes and avoid overlaps/holes
            # due to rounding errors
            x1 = xmin + (col * hOverlay)    # far left
            x2 = x1 + (xVertexHi - xVertexLo)  # left
            x3 = xmin + (col * hOverlay) + hSpacing  # right
            x4 = x3 + (xVertexHi - xVertexLo)  # far right
    
            for row in range(rows):
                if (col % 2) == 0:
                    y1 = ymax + (row * vOverlay) - (((row * 2) + 0) * halfVSpacing)  # hi
                    y2 = ymax + (row * vOverlay) - (((row * 2) + 1) * halfVSpacing)  # mid
                    y3 = ymax + (row * vOverlay) - (((row * 2) + 2) * halfVSpacing)  # lo
                else:
                    y1 = ymax + (row * vOverlay) - (((row * 2) + 1) * halfVSpacing)  # hi
                    y2 = ymax + (row * vOverlay) - (((row * 2) + 2) * halfVSpacing)  # mid
                    y3 = ymax + (row * vOverlay) - (((row * 2) + 3) * halfVSpacing)  # lo
    
                geoms.append((
                    (x1, y2),
                    (x2, y1), (x3, y1), (x4, y2), (x3, y3), (x2, y3),
                    (x1, y2)
                ))
                
        if self._cut == True:
            hexData = geopandas.GeoDataFrame(
                index=[i for i in range(len(geoms))],
                geometry=pandas.Series(geoms).apply(lambda x: shapely.geometry.Polygon(x)),
                crs=gdf.crs)
            return geopandas.overlay(gdf, hexData, how='intersection')
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
        cls : TYPE
            DESCRIPTION.
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
        sanitized : TYPE
            DESCRIPTION. A collection of polygon geometry as a geodataframe object

        """
        cls._gdf = gdf
        cls._tolerance = tolerance
        if cls._tolerance < 11:
            raise ValueError(f" -- Tolerance can not be less than 11 (1.1 centimeter). \
                             Please assign tolerance value more than 11")
        
        if cls._gdf.crs == 'EPSG:4326':
        
            geoms = [i.buffer(0.00000001).buffer(-cls._tolerance/1e9)
                     for i in [items for sublist in [list(i)
                                                     for i in cls._gdf.geometry if i.geom_type == 'MultiPolygon']
                               for items in sublist]
                     + [i for i in cls._gdf.geometry if i.geom_type == 'Polygon']]
        
            sanitized = geopandas.GeoDataFrame(crs = 'EPSG:4326', geometry = geoms)
        else:
            toGdf = cls._gdf.to_crs('EPSG:4326')
            geoms = [i.buffer(0.00000001).buffer(-cls._tolerance/1e9)
                     for i in [items for sublist in [list(i)
                                                     for i in toGdf.geometry if i.geom_type == 'MultiPolygon']
                               for items in sublist]
                     + [i for i in toGdf.geometry if i.geom_type == 'Polygon']]    
        
            sanitized = geopandas.GeoDataFrame(crs = 'EPSG:4326', geometry = geoms)
        
        return sanitized