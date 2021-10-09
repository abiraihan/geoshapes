# geoshapes
Geospatial Experimental Design with GIS Analytics

### Install pypi library
```python
pip install geoshapes
```
#### Example
```python
import string, shapely, geosolution, geopandas
point = shapely.geometry.Point(0,0)
   
c = geosolution.splitShape.splitLatin(point, 25)
ft = geopandas.GeoDataFrame(geometry = c, crs = 'EPSG:4326')

ft['ids'] = range(len(ft))
ft['Group']= ft.apply(lambda row : string.ascii_uppercase[int(row.ids)], axis = 1)
ax = ft.plot(figsize=(15, 10), alpha=0.4, edgecolor='k')
ft.plot(column='Group', ax=ax, linewidth=9, cmap='tab20');
ft.apply(lambda x: ax.annotate(s=f"{x.Group}",
                               xy=x.geometry.centroid.coords[0],
                               weight='bold',
                               ha='center',
                               va='center',
                               size=20),axis=1)
```
### Split a circle geometry with defined indentifier as a treatment plot
<p align="center">
<img src="https://github.com/abiraihan/geoshapes/blob/master/images/splitCircle.png" width="600">
</p>

License
----
MIT
