**gridShape**
==============
gridShape module functionality involved to create grid geometry for a given boundary geometry.

squareGrid
------------

*Create square grid polygon geometry for given boundary geometry*

:module: geoshapes.gridShape.squareGrid

.. function:: squareGrid(boundary, length:int, cut:bool = False)

   :param boundary: Boundary shapefile or geopandas GeoDataFrame type
   :type boundary: str / geopandas.GeoDataFrame
   :param length: length of the side of a square grid in meters
   :type length: int
   :param cut: intersects grid polygon with boundary area
   :type cut: True
   :return: Square grid of the boundary
   :rtype: geopandas.GeoDataFame
    
.. container:: header

    **Code Block**

.. code-block:: python

  import geoshapes, gropandas
  
  fileData = geopandas.read_file("./polygonShapefileData.shp")
  bounds = fileData.geometry[0]
  dr = geopandas.GeoDataFrame(crs = 'EPSG:3857', geometry = [bounds])
  dr = dr.to_crs('EPSG:4326')
  fr = geoshapes.gridShape.squareGrid(dr, 10, cut = True)
  fr.plot(cmap = 'tab20')
  
.. container:: header

        *Output Map*
.. image:: ../images/squareGrid.png
   :scale: 80 %
   :alt: squareGrid Output
   :align: center

----------------------------------------------------------------------------------------------------