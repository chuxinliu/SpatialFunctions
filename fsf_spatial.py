'''
Functions to get spatial information:
1. (parcel level) distance to coastline
2. (parcel level) elevation 

Function to merge BK and FSF using nearest point;
'''


# import libraries
import os
from math import *
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge 

import seaborn as sns
import matplotlib.pyplot as plt

from haversine import haversine, Unit
from shapely import geometry
from shapely.ops import nearest_points
from shapely import wkt
from shapely.geometry import Point, Polygon
from shapely.wkt import loads
import contextily as cx
from scipy.spatial import cKDTree

# ------------------------------------------------------------------------

# get coastline for a state based on it's max/min lon/lat
def get_state_coastline(df, coastline, lat_column, lon_column, margin):
	'''
	Input:
	df: pandas df that has lon and lat
	coastline: geopandas object of the world's coastline
	lat_column, lon_column: string, names of dataframe columns
	margin: float, extra degrees to cover the state
	----------------------------------------------------
	Output: geo pandas of coastline that covers this state
	'''
	
	lo1,lo2 = min(df[lon_column]),max(df[lon_column])	
	la1,la2 = min(df[lat_column]),max(df[lat_column])

	p1 = geometry.Point(lo1-margin,la1-margin)
	p2 = geometry.Point(lo1-margin,la2+margin)
	p3 = geometry.Point(lo2+margin,la2+margin)
	p4 = geometry.Point(lo2+margin,la1-margin)

	points = [p1,p2,p3,p4,p1]

	shapely_poly = geometry.Polygon([[p.x, p.y] for p in points])
	pd_poly = pd.DataFrame( data = {'aoi':[1], 
									'wkt':[shapely_poly.wkt]} )
	gpd_poly = gpd.GeoDataFrame(pd_poly, 
								crs=4326, 
								geometry=[loads(pgon) for pgon in pd_poly.wkt])
	coastline_aoi = gpd.sjoin(coastline, gpd_poly, how = 'inner')

	return coastline_aoi


# get a tuple of (closest point, haversine distance) between a point and a column of string lines
def point_to_coast_min_dist(point, geodf, unit):
	'''
	Input: 
	point: one point's geo
	geodf: geo dataframe which geometry column is a bunch of coastline stringlines
	unit: meters or miles
	-----------------------------------------------------------------------------
	Output: (a tuple) 
	- two points as a list
	- smallest distance in unit of meters as a float
	'''
	closest = []
	min_dist = 10**9
	
	for i in range(geodf.count()[0]):
		two_points = nearest_points(point, geodf.iloc[i]['geometry'])
		dist = two_points[0].distance(two_points[1])
	
	if dist<min_dist:
	  	closest = [two_points[0], two_points[1]]
	  	min_dist = dist

	if closest:
		return closest, haversine((closest[0].y, closest[0].x), (closest[1].y, closest[1].x), unit=unit)
	else: 
		print("no closest distance for this point", point)
		return None

# get the smallest haversine distance to coastline
def points_to_coastlines(gdf1, gdf2, unit):
	'''
	input: two geo dataframe, gdf1 has all points, gdf2 has stringlines
	unit: Unit.MILES or Unit.METERS
	---------
	!!! inplace change, no return !!!
	add an new column to gdf1 called "dist_to_coast"
	'''
	closest_dist = []
	for i in range(gdf1.count()[0]):
		point = gdf1.iloc[i]['geometry']
		min_dist = point_to_coast_min_dist(point, gdf2, unit)[1]
		closest_dist.append(min_dist)

	array_dist = np.array(closest_dist)
	gdf1['dist_to_coast'] = array_dist.tolist()

# merge DEM rasterio objects
def merge_dem_rio_into_1rio(list_rio):
	'''
	input: list_rio: a list of raster io objects
	output: 1 rio 
	'''
	merged, transf = merge(list_rio)

	output_meta = n38.meta.copy()
	output_meta.update( {"driver": "GTiff",
						"height": merged.shape[1],
	    				"width": merged.shape[2],
	    				"transform": transf} )
	with rasterio.open('merged.tif', "w", **output_meta) as m:
		m.write(merged)

	merged_rio = rasterio.open('merged.tif')
	return merged_rio



# get elevation using rasterio, DEM (as in tif file)
def get_elevation(df, lat_column, lon_column, dem):
	"""
	Input:
	dem: rasterio object 
	df:  pandas or geopandas dataframe
	lat_column, lon_column: string, names of dataframe columns
	----------------------------------------------------------
	!!! inplace change, no return !!!
	df['elev_meters'] column is added to df
	"""

	elev_meters = []
	for index, row in df.iterrows():
		# extra x and y
		x, y = row[lon_column], row[lat_column]

		# take elevation value from dem file
		row, col = dem.index(x, y)
		elev = dem.read(1)[row, col]
		elev_meters.append(elev)

	df['elev_meters'] = elev_meters


def ckdnearest(gdfA, gdfB):
	'''
	Input:
	gdfA: the dataframe to keep
	gdfB: the dataframe to find nearest point's value
	
	Output:
	gpd: geodataframe merged with nearest point and its distance
	'''

	nA = np.array(list(gdfA.geometry.apply(lambda x: (x.x, x.y))))
	nB = np.array(list(gdfB.geometry.apply(lambda x: (x.x, x.y))))

	btree = cKDTree(nB)
	dist, idx = btree.query(nA, k=1)
	gdfB_nearest = gdfB.iloc[idx].drop(columns="geometry").reset_index(drop=True)

	gdf = pd.concat([
		gdfA.reset_index(drop=True),
		gdfB_nearest,
		pd.Series(dist, name = 'dist_to_nearest')
		]
		,axis=1)

	return gdf 

#########################################################################
#########################################################################

# read bk
bk = pd.read_csv('attom_geocode_10_DE.csv')
bk = bk[['sa_property_id','x','y']]
geobk = gpd.GeoDataFrame(bk, geometry=gpd.points_from_xy(bk.x, bk.y), crs='EPSG:4326')

# read FSF data
de = pd.read_stata('FSF_10_info.dta')
de = de[['fsid','lon','lat','floodfactor','avm']]
geode = gpd.GeoDataFrame(de, geometry=gpd.points_from_xy(de.lon, de.lat))

# read coastline data
coastline = gpd.read_file('tl_2019_us_coastline.shp').to_crs(epsg = 4326)

# read elevation raster data
n38 = rasterio.open('n38_w076_1arc_v3.tif', crs='EPSG:4326')
n39 = rasterio.open('n39_w076_1arc_v3.tif', crs='EPSG:4326')

# ------------------------------------------------------------------------

# get coastline's "area of interest" geodataframe
coastline_aoi = get_state_coastline(de, coastline, 'lat', 'lon', 0.5)

# (inplace) get the smallest haversine distance to coastline
points_to_coastlines(geode, coastline_aoi, Unit.METERS)

# merge dem rasterio objects into 1
list_rio = [n38,n39]
de_rio = merge_dem_rio_into_1rio(list_rio)

# (inplace) get elevation
get_elevation(geode, 'lat', 'lon', de_rio)

geode.to_csv('geode.csv')

print("FSF‘s data frame running done")

# ------------------------------------------------------------------------

bk_nearest_fsf = ckdnearest(geobk, geode)

bk_nearest_fsf.to_csv('bk_nearest_fsf.csv')

print("Merge BK+FSF running done")






