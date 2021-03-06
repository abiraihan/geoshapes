**checkShape**
==============
checkShape module functionality involved to check geometry validity, overlaps and geometry type.

checkPolygon
------------

*Chcek polygon geometry validity and returns only valid polygon geometry*

:module: geoshapes.checkShape.checkPolygon

.. function:: checkPolygon(data:str)

   :param data: data file or directory path
   :type data: str
   :return: geopandas GeoDataFame
   :rtype: geopandas.GeoDataFame
    
.. container:: header

    **Code Block**

.. code-block:: python

  import geoshapes
  
  path = "./polygonShapefileData"
  dataFiles = geoshapes.checkShape.checkPolygon(path)
  
.. container:: header

        *Output*
   -- Found No Issue : Performed Topology Check for the geometry


----------------------------------------------------------------------------------------------------

checkOverlaps
------------

*To identify and check the topology of polygon data layer*

:module: geoshapes.checkShape.checkOverlaps

.. function:: checkOverlaps(geomPaths:str, pathType:str, tolerance:float)

   :param geomPaths: shapefile or directory path
   :type geomPaths: str
   :param pathType: default 'directoryPath'
   :type pathType: str
   :param tolerance: default is 0.005
   :type tolerance: float, optional
   :return: geopandas GeoDataFame
   :rtype: geopandas.GeoDataFame
    
.. container:: header

    **Code Block**

.. code-block:: python

  import geoshapes
  
  path = "./polygonShapefileData"
  dataFiles = geoshapes.checkShape.checkOverlaps(path)
  
.. container:: header

        *Output*
    -- Please check the validity (Topology) of geometry data layer
    -- Found overlapped geometry into the data layer

----------------------------------------------------------------------------------------------------

checkGeometry
-------------

*To identify and check the topology of polygon data layer*

:module: geoshapes.checkShape.checkGeometry

.. function:: checkGeometry(fileObjects:geopandas.GeoDataFrame, geomType:str)

   :param fileObjects: geopandas GeoDataFrame File Type
   :type geomPaths: geopandas.GeoDataFrame
   :param geomType: geometry type i.e. Point, Polygon, LineString
   :type pathType: str
   :return: boolean
   :rtype: boolean
    
.. container:: header

    **Code Block**

.. code-block:: python

  import geoshapes
  
  path = "./polygonShapefileData"
  dataFiles = geoshapes.checkShape.checkPolygon(path)
  isPolygon = geoshapes.checkShape.checkGeometry(dataFiles, 'Polygon')
  print(isPolygon)
  
.. container:: header

        *Output*
    True

----------------------------------------------------------------------------------------------------