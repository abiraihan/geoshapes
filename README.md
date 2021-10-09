# geoshapes
Geospatial Experimental Design with GIS Analytics

### Install pypi library
```python
pip install geoshapes
```
#### Example
```python
 1  import string, shapely, geoshapes, geopandas
 2  pointLocation = shapely.geometry.Point(0,0)
 3  polygonList = geoshapes.splitShape.splitCircle(geoms = pointLocation,
 4                                                   circleRadius = 500,
 5                                                   incrementdegree = 45,
 6                                                   clipInterior = True,
 7                                                   innerWidth = 100,
 8                                                   getGeom = 'Both'
 9                                                   )
10  gdf = geopandas.GeoDataFrame(geometry = polygonList, crs = 'EPSG:3857')
11  gdf['ids'] = range(len(gdf))
12  gdf['Group']= gdf.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
13  ax = gdf.plot(figsize=(15, 10), alpha=0.0, edgecolor='k')
14  gdf.plot(column='Group', ax=ax, linewidth=9, cmap='tab20');
15  gdf.apply(lambda x: ax.annotate(s=f"Group : {x.Group}{x.ids}",
16                                  xy=x.geometry.centroid.coords[0],
17                                  weight='bold', ha='center',
18                                  va='center', size=10),axis=1
19                                  )
```
### Split a circle geometry with defined indentifier as a treatment plot
<p align="center">
<img src="https://github.com/abiraihan/geoshapes/blob/master/images/splitCircle.png" width="600">
</p>

License
----
MIT
