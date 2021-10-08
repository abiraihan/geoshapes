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
c = splitShape.splitCircle(pointShape, 150, 45, clipInterior = False, innerWidth = 50, getGeom = 'Both')
ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:3857')
ft['ids'] = range(len(ft))
ft['Group']= ft.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
ft.plot(cmap = 'tab20')
```
### Split a circle geometry with defined indentifier as a treatment plot
<p align="center">
<img src="https://github.com/abiraihan/geoshape/blob/master/images/splitCircle.png" width="600">
</p>

License
----
MIT
