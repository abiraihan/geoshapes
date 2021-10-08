# geoshape
Geospatial Experimental Design with GIS Analytics

### Install pypi library
```python
pip install geoshape
```
#### Example
```python
import shapely
import geoshape
pointShape = shapely.geometry.Point(0.0, 0.0)
pointPolys = pointShape.buffer(300)
circleGeoms = geoshape.splitCircle(pointPolys)
```

License
----
MIT
