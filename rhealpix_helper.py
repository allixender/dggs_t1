from rhealpixdggs.dggs import *
import pandas as pd
import geopandas as gpd
from pyproj import Transformer
from shapely.ops import transform
try:
    import rasterio
except ImportError:
    print("rasterio not available")

import time
import sys
import os
import matplotlib.pyplot as plt
sys.path.append('..')
from shapely.geometry import Polygon, Point, box
from pandas.core.common import flatten
from math import radians, sin, cos, asin, sqrt

def __haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r * 1000


def __lonlat_to_latlon(lonlat_array):
    latlon_array = []
    for vertex in lonlat_array:
        latlon_array.append((vertex[1],vertex[0]))
    return latlon_array

def __cell_to_geometry(cell):
    geom = None
    try:
        # geom =  Polygon(__lonlat_to_latlon(cell.boundary(n=2,plane=False)))
        # gdf['geometry'] = gdf['cell_id'].apply(lambda x: Polygon(x.boundary(n=2,plane=False)))
        geom = Polygon(cell.boundary(n=2,plane=False))
    except:
        print(f'internal rhealpix error with cell.boundary method for {str(cell)}')
    return geom

def create_rhpix_geometry(df):

    gdf = gpd.GeoDataFrame(df.copy())
    gdf['geometry'] = gdf['cell_id'].apply(__cell_to_geometry)

    gdf.crs = 'EPSG:4326'
    gdf['cell_id'] = gdf['cell_id'].apply(lambda x: str(x))

    return gdf


def get_rhpix_cells(res, extent=None):
    rdggs = WGS84_003
    if extent:
        se = (extent[1], extent[2])
        nw = (extent[3], extent[0])
        set_hex = list(flatten(rdggs.cells_from_region(res, se, nw, plane=False)))
    else:
        set_hex = [x for x in rdggs.grid(res)]

    df = pd.DataFrame({"cell_id": set_hex})

    return df


def create_rhpix_geom_cells_global(resolutions, table, export_type, db_engine=''):
    """Create geometry for rhpix cells globally for given resolutions

        Parameters:
        db_engine (sqlalchemy.engine): sqlalchemy database engine
        resolutions(array): array of integer h3 resolution levels
        table(string): table name for postgres database
        export_type(string): where to export 'geojson' or 'postgres'

        Returns:
        none
    """
    rdggs = WGS84_003
    transformer = Transformer.from_crs("epsg:4326", 'proj=rhealpix')
    for res in resolutions:

        gdf = gpd.GeoDataFrame({'cell_id':[x for x in rdggs.grid(res)]})
        gdf['geometry'] = gdf['cell_id'].apply(lambda x: Polygon(x.boundary(n=10,plane=False)))
        gdf.crs = 'EPSG:4326'
        gdf['cell_id'] = gdf['cell_id'].apply(lambda x: str(x))
        gdf['area'] = gdf['geometry'].apply(lambda x: transform(transformer.transform, x).area)

        print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))

        if export_type == 'postgres':

            gdf.to_postgis(table + str(res), db_engine, if_exists='replace')
            print('finish import to db {} {}'.format(res, time.asctime(time.localtime(time.time()))))


        elif export_type == 'geojson':

            gdf.to_file("{}{}.geojson".format(table, res), driver='GeoJSON')
            print('finish import to geojson {} {}'.format(res, time.asctime(time.localtime(time.time()))))


def raster_to_rhpix(raster_path, value_name, cell_min_res, cell_max_res, extent=None, pix_size_factor=3):
    """Load raster values into h3 dggs cells

    Parameters:
    raster (string): path to raster file for uploading
    resolutions (string): srs epsg code of raster's territory coordinate system
    table (string): name of a value to be uploaded into dggs cells
    cell_min_res (integer): min h3 resolution to look for based on raster cell size
    cell_max_res (integer): max h3 resolution to look for based on raster cell size
    extent (list): Extent as array of 2 lon lat pairs to get raster values for
    pix_size_factor (pinteger): how times smaller h3 hex size should be comparing with raster cell size

    Returns:
    Pandas dataframe
   """
    # Open raster
    rs = rasterio.open(raster_path)
    rdggs = WGS84_003
    # Get extent to fill with rhealpix cells

    if extent:
        nw = (extent[1], extent[2])
        se = (extent[3], extent[0])
    else:
        extent = gpd.GeoSeries(box(rs.bounds.left, rs.bounds.bottom, rs.bounds.right, rs.bounds.top)).__geo_interface__

    # Get resolution:edge lenght in m dict
    resolutions = {}
    for i in range(cell_min_res, cell_max_res, 1):
        resolutions[i] = rdggs.cell_width(i)

    # Get two points on borders of neighbour pixels in raster
    x1 = rs.transform[2]
    y1 = rs.transform[5]
    x2 = rs.transform[2] + rs.transform[0]
    y2 = rs.transform[5] - rs.transform[4]

    # Get pixel size from projected src
    #     transformer = Transformer.from_crs("epsg:4326", 'proj=isea')
    #     size = Point(transformer.transform(y1, x1)).distance(Point(transformer.transform(y2, x1)))
    #     print(f'Projected size{size}')

    # Get pixel size from haversine formula

    size = __haversine(x1, y1, x1, y2)

    print(f"Raster pixel size {size}")

    # Get raster band as np array
    raster_band_array = rs.read(1)

    # Get h3 resolution for raster pixel size
    for key, value in resolutions.items():
        if value < size / pix_size_factor:
            resolution = key
            break
    print(resolution)

    # Create dataframe with cell_ids from extent with given resolution
    print(f"Start filling raster extent with rhealpix indexes at resolution {resolution}")
    df = pd.DataFrame({'cell_id': list(flatten(rdggs.cells_from_region(resolution, nw, se, plane=False)))})

    # Get raster values for each cell_id
    print(f"Start getting raster values for cells at resolution {resolution}")
    df[value_name] = df['cell_id'].apply(lambda x: raster_band_array[rs.index(x.centroid(plane=False)[1], x.centroid(plane=False)[0])])

    # Drop nodata
    df = df[df[value_name] != rs.nodata]

    return df
