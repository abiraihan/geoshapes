# geoshape
Geospatial Experimental Design with GIS Analytics

### Install geoshape library
```
pip install geoshape
```
#### Example
```
import shapely
import geoshape
pointShape = shapely.geometry.Point(0.0, 0.0)
pointPolys = pointShape.buffer(300)
circleGeoms = geoshape.splitCircle(pointPolys)
```
