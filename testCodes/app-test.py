import shapely, geopandas
from geoshapes import splitShape

pointLocation = shapely.geometry.Point(0.0,0.0)
circle = splitShape.splitCircle(
    geoms = pointLocation,
    circleRadius = 50,
    incrementDegree = 5,
    clipInterior = True,
    innerWidth = 30,
    getGeom = 'Both'
    )

gdf = geopandas.GeoDataFrame(geometry = circle)
print(f"Total row number is geodataframe is {len(gdf)}")

areas = sum([i.area*1e10/4046.86 for i in gdf.geometry])
print(f"Total ares of the geometry is {round(areas, 2)} acre")