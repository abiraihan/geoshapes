.. image:: https://www.repostatus.org/badges/latest/active.svg
   :alt: Project Status: Active – The project has reached a stable, usable state and is being actively developed.
   :target: https://www.repostatus.org/#active

.. image:: https://badge.fury.io/py/geoshapes.svg
    :target: https://badge.fury.io/py/geoshapes

**mergeShape**
==============
mergeShape module functionality involved to merge splited geometry for a defined number of geometry.

mergePolygon
------------

*Merge splited polygon geometry for a defined number of geometry*

:module: geoshapes.mergeShape.mergePolygon

.. function:: mergePolygon(geomData, mergedPoly:int)

   :param geomData: list of shapely polygon geometry or geopandas GeoDataFrame
   :type geomData: list/geopandas.GeoDataFrame
   :param mergedPoly: least Number of Merged geometry.
   :type mergedPoly: int
   :return: A geopandas GeoDataFrame
   :rtype: geopandas.GeoDataFame
    
.. container:: header

    **Code Block**

.. code-block:: python

  import geoshapes, gropandas
  
  fileData = geopandas.read_file("./polygonShapefileData.shp")
  bounds = fileData.geometry[0]
  dr = geopandas.GeoDataFrame(crs = 'EPSG:3857', geometry = [bounds])
  dr = dr.to_crs('EPSG:4326')
  square = geoshapes.gridShape.squareGrid(dr, 10, cut = True)
  square.plot(
        figsize=(8, 7),
        alpha=0.3,
        edgecolor='k'
        )
  
.. container:: header

        *Output Map*
.. image:: ../docs/images/squareGrid.png
   :scale: 80 %
   :alt: squareGrid Output
   :align: center

----------------------------------------------------------------------------------------------------