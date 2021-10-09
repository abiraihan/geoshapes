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

  import shapely, geoshapes, geopandas
  
  path = "./polygonShapefileData"
  dataFiles = geoshapes.checkShape.checkPolygon(path)
  
.. container:: header

        *Output*
.. note::
   -- Found No Issue : Performed Topology Check for the geometry


----------------------------------------------------------------------------------------------------