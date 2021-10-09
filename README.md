# geoshapes
Geospatial Experimental Design with GIS Analytics

### Install library
```python
pip install geoshapes
```
#### Example
#####        geoshapes.splitShape.splitCircle
######         - Split a circle geometry with defined indentifier as a treatment plot
```python
import string, shapely, geoshapes, geopandas
pointLocation = shapely.geometry.Point(0,0)
polygonList = geoshapes.splitShape.splitCircle(
    geoms = pointLocation,
    circleRadius = 500,
    incrementDegree = 45,
    clipInterior = True,
    innerWidth = 100,
    getGeom = 'Both'
    )
gdf = geopandas.GeoDataFrame(
    geometry = polygonList,
    crs = 'EPSG:3857'
    )
gdf['ids'] = range(len(gdf))
gdf['Group']= gdf.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
ax = gdf.plot(figsize=(15, 10), alpha=0.0, edgecolor='k')
gdf.plot(column='Group',
         ax=ax, linewidth=9,
         cmap='tab20');
gdf.apply(lambda x: ax.annotate(
    s=f"Group : {x.Group}{x.ids}",
    xy=x.geometry.centroid.coords[0],
    weight='bold', ha='center',
    va='center', size=10),axis=1
    )
```

  Figure : splitsCircle
_____

<p align="center">
<img src="https://github.com/abiraihan/geoshapes/blob/master/images/splitCircle.png" width="600">
</p>

#####        geoshapes.splitShape.splitLatin
######         - Split a square geometry with defined indentifier as a latin square treatment plot
```python
import string, shapely, geoshapes, geopandas
point = shapely.geometry.Point(0,0)
latinShpae = geoshapes.splitShape.splitLatin(point, 25)
featureGeoms = geopandas.GeoDataFrame(
    geometry = latinShpae,
    crs = 'EPSG:4326'
    )
featureGeoms['ids'] = range(len(featureGeoms))
featureGeoms['Group']= featureGeoms.apply(
    lambda row : string.ascii_uppercase[int(row.ids)],
    axis = 1
    )
ax = featureGeoms.plot(
    figsize=(15, 10),
    alpha=0.4, edgecolor='k')
featureGeoms.plot(
    column='Group', ax=ax,
    linewidth=9, cmap='tab20'
    );
featureGeoms.apply(
    lambda x: ax.annotate(
        s=f"{x.Group}",
        xy=x.geometry.centroid.coords[0],
        ha='center',
        va='center',
        size=20),axis=1
    )
```

  Figure : splitLatin
_____

<p align="center">
<img src="https://github.com/abiraihan/geoshapes/blob/master/images/latinSquare.png" width="600">
</p>

License
----
MIT
