.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/abiraihan/geoshapes/b8336cb953c2060a7b0209b7732db6859e26ea80?urlpath=lab%2Ftree%2Fexample%2FUntitled.ipynb

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

    import string, shapely, geoshapes, geopandas
    polys = shapely.geometry.Polygon([(0, 0), (0,5), (16, 5), (16, 0)])
    splitGeometry = geoshapes.splitShape.splitGeom(polys, 50, rotation = 120)
    
    #Input Polygon geometry for mergePolygon function
    splitGeometry['ids'] = range(len(splitGeometry))
    ax1 = splitGeometry.plot(figsize = (7,5), alpha = 0.9, cmap = 'Spectral', edgecolor = 'k', linewidth = 1)
    ax1.set_title('splitShape.splitGeom | Splited Polygon : 50')
    splitGeometry.apply(
        lambda x: ax1.annotate(
            s=f"{x.ids+1}",
            xy=x.geometry.centroid.coords[0],
            ha='center',
            va='center',
            size=5
            ),
        axis=1
        )
    
    #Output Polygon geometry for mergePolygon function
    geoData = geoshapes.mergeShape.mergePolygon(splitGeometry, 3)
    geoData['ids'] = range(len(geoData))
    ax = geoData.plot(figsize = (7,5), alpha = 0.9, cmap = 'Spectral', edgecolor = 'k', linewidth = 2)
    ax.set_title('mergeShape.mergePolygon | Merged Polygon : 3')
    geoData.apply(
        lambda x: ax.annotate(
            s=f"{x.ids+1}",
            xy=x.geometry.centroid.coords[0],
            ha='center',
            va='center',
            size=10
            ),
        axis=1
        )

.. container:: header

        *Input Polygon*
        
.. image:: ../docs/images/mergePolygonInput.png
   :scale: 80 %
   :alt: splitedPolygon Input
   :align: center
   
.. container:: header

        *Output Map*

.. image:: ../docs/images/mergePolygon.png
   :scale: 80 %
   :alt: mergePolygon Output
   :align: center

----------------------------------------------------------------------------------------------------