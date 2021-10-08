# geoshape
Geospatial Experimental Design with GIS Analytics

### Install pypi library
```python
pip install geoshape
```
#### Example
```python
import string, shapely, geoshape, geopandas
pointShape = shapely.geometry.Point(0.0, 0.0)
c = splitShape.splitCircle(pointShape, 3188, 5, clipInterior = True, innerWidth = 2788, getGeom = 'Outer')
ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:3857')
ft['ids'] = range(len(ft))
ft['Group']= ft.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
ft.plot(cmap = 'tab20')
```

License
----
MIT
