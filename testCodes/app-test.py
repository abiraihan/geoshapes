import shapely, geopandas
from geoshapes import splitShape

pointLocation = shapely.geometry.Point(0.0,0.0)
circle = splitShape.splitCircle(
    geoms = pointLocation,
    circleRadius = 50,
    incrementDegree = 5,
    clipInterior = True,
    innerWidth = 30,
    getGeom = 'Both')
gdf = geopandas.GeoDataFrame(geometry = circle)

print(f"Total row number is geodataframe is {len(gdf)}")

areas = sum([i.area*1e10/4046.86 for i in gdf.geometry])
print(f"Total ares of the geometry is {round(areas, 2)} acre")

fileDatas = geopandas.GeoDataFrame()
for i in gdf.geometry:
    f = splitShape.splitGeom(i, 5, rotaion = 100)
    ft = geopandas.GeoDataFrame(geometry = f, crs = 'EPSG:4326')
    fileDatas  =fileDatas.append(ft)

fileDatas.reset_index(drop = True, inplace = True)
print(f"Total row number is geodataframe is {len(fileDatas)}")