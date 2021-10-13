# geoshapes
[![PyPI version](https://badge.fury.io/py/geoshapes.svg)](https://badge.fury.io/py/geoshapes)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/abiraihan/geoshapes.git/master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5559438.svg)](https://doi.org/10.5281/zenodo.5559438)

Geospatial Experimental Design with GIS Analytics.
To split Polygon geometry into different shape that required to create
experimental plot / Trial design. It creates/re-creates experimental
plot design to summerize data to geospatially enabled space for next-step
mathematical modelling.


### Run example code on mybinder.org
You can run geoshapes example code on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/abiraihan/geoshapes/8d441eef49cd387980a86ec230a84fde012390a3?urlpath=lab%2Ftree%2Fexample%2FsplitShape.ipynb) without installing any pthon library in your python environment.
It doesn't requie any dependency to install and you can run all those example from the mybinder.org


#### * *<a href="./docs/usage.rst">Required Library</a> | Required python library prior to install geoshapes*

### Install geoshapes
```python
pip install geoshapes
```

### Module Summery:

#### *1. <a href="./docs/splitShape.rst">geoshapes.splitShape</a> | Splits/create geometry for experimetal design.*
#### *2. <a href="./docs/checkShape.rst">geoshapes.checkShape</a> | Check geometry validity & return/report valid geometry.*
#### *3. <a href="./docs/gridShape.rst">geoshapes.gridShape</a> | Create grid for a given boundary geometry.*


### Example
##### geoshapes.splitShape.splitCircle
###### - Split a circle geometry with defined indentifier as a treatment plot
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
gdf['Group']= gdf.apply(
    lambda row : string.ascii_uppercase[int(row.ids)],
    axis = 1)
ax = gdf.plot(figsize=(15, 10),
alpha=0.0, edgecolor='k')
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

![Split CIrcle](https://github.com/abiraihan/geoshapes/blob/master/docs/images/splitCircle.png)
_____



##### geoshapes.splitShape.splitLatin
###### - Split a square geometry with defined indentifier as a latin square treatment plot
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

![Latin Square](https://github.com/abiraihan/geoshapes/blob/master/docs/images/latinSquare.png)
_____


License
----
MIT
