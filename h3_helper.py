from h3 import h3
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point, box
from shapely.ops import transform
import time
from pyproj import Transformer
import rasterio
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


def create_h3_geometry(df):
    gdf = gpd.GeoDataFrame(df)
    gdf['geometry'] = df['cell_id'].apply(lambda x: Polygon(h3.h3_to_geo_boundary(x,geo_json=True)))
    gdf.crs = 'EPSG:4326'
    return gdf


def create_h3_geom_cells(extent, resolutions, table, export_type, db_engine):
    """Create geometry for h3 cells in given extent for given resolutions levels

    Parameters:
    extent (array): array of lat lon extent pairs for covering with h3 cells
    resolutions(array): array of integer h3 resolution levels
    table(string): table name for postgres database
    db_engine (sqlalchemy.engine): sqlalchemy database engine
    export_type(string): where to export 'geojson' or 'postgres'

    Returns:
    none
   """
    extent = gpd.GeoSeries(box(extent[0], extent[1], extent[2], extent[3])).__geo_interface__
    for res in resolutions:
        calc_time = time.time()
        print(f'start caclulating resolution {res}')
        set_hex = list(h3.polyfill(extent['features'][0]["geometry"], res=res))

        print(f'finish caclulating resolution {res} in {str(round(time.time() - calc_time, 2))} seconds')

        if export_type == 'postgres':
            geom_time = time.time()
            gdf = gpd.GeoDataFrame({"cell_id": set_hex})
            gdf['geometry'] = gdf["cell_id"].apply(lambda x: (Polygon(h3.h3_to_geo_boundary(x,geo_json=True))))

            print(f'finish caclulating geometry fo res {res} in {str(round(time.time() - geom_time, 2))} seconds')

            import_time = time.time()
            gdf.to_postgis(table + str(res), db_engine, if_exists='replace')

            print(f'finish import to db {res} in {str(round(time.time() - import_time, 2))} seconds')

        elif export_type == 'geojson':
            transformer = Transformer.from_crs("epsg:4326", 'proj=isea')
            gdf = gpd.GeoDataFrame({"cell_id": set_hex})

            gdf['geometry'] = gdf["cell_id"].apply(lambda x: Polygon(h3.h3_to_geo_boundary(x, True)))
            gdf['area'] = gdf["geometry"].apply(lambda x: transform(transformer.transform, x).area)
            gdf.to_file("{}{}.geojson".format(table, res), driver='GeoJSON')
            print('finish import to geojson {} {}'.format(res, time.asctime(time.localtime(time.time()))))


def create_h3_geom_cells_global(resolutions, table, export_type, db_engine=''):
    """Create geometry for h3 cells globally for given resolutions

        Parameters:
        db_engine (sqlalchemy.engine): sqlalchemy database engine
        resolutions(array): array of integer h3 resolution levels
        table(string): table name for postgres database
        export_type(string): where to export 'geojson' or 'postgres'

        Returns:
        none
    """
    for res in resolutions:
        set_hex_0 = list(h3.get_res0_indexes())
        set_hex = []
        if res == 0:
            set_hex = set_hex_0
        else:
            for i in set_hex_0:
                set_hex.extend(list(h3.h3_to_children(i, res)))
        if export_type == 'postgres':
            gdf = pd.GeoDataFrame({"cell_id": set_hex})
            gdf['geometry'] = gdf["cell_id"].apply(lambda x:(Polygon(h3.h3_to_geo_boundary(x, geo_json=True)).wkb))

            print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))

            gdf.to_postgis(table + str(res), db_engine, if_exists='replace')
            print('finish import to db {} {}'.format(res, time.asctime(time.localtime(time.time()))))

        elif export_type == 'geojson':
            transformer = Transformer.from_crs("epsg:4326", 'proj=isea')
            gdf = gpd.GeoDataFrame({"cell_id": set_hex})
            gdf['geometry'] = gdf.cell_id.apply(lambda x: Polygon(h3.h3_to_geo_boundary(x, geo_json=True)))
            print('finish caclulating geometry {} {}'.format(res, time.asctime(time.localtime(time.time()))))
            gdf['area'] = gdf.geometry.apply(lambda x: transform(transformer.transform, x).area)
            gdf.to_file("{}{}.geojson".format(table, res), driver='GeoJSON')
            print('finish import to geojson {} {}'.format(res, time.asctime(time.localtime(time.time()))))


def get_h3_cells(res, extent=None):
    
    """Get h3 cells for given resolution

    Parameters:
    res (int): h3 resolution 
    extent (list): Extent as array of 2 lon lat pairs to get raster values for
    Returns:
    Pandas dataframe
   """
    if extent:
        set_hex = list(h3.polyfill_geojson(extent, res=res))
    else:    
        set_hex_0 = list(h3.get_res0_indexes())
        set_hex = []
        if res == 0:
            set_hex = set_hex_0
        else:
            for i in set_hex_0:
                set_hex.extend(list(h3.h3_to_children(i, res)))
    df = pd.DataFrame({"cell_id": set_hex})
    
    return df


def raster_to_h3(raster_path, value_name, cell_min_res, cell_max_res, extent=None, pix_size_factor=3):
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

    # Get extent to fill with h3 hexes
    if extent:
        extent = gpd.GeoSeries(box(extent[0], extent[1], extent[2], extent[3])).__geo_interface__
    else:
        extent = gpd.GeoSeries(box(rs.bounds.left, rs.bounds.bottom, rs.bounds.right, rs.bounds.top)).__geo_interface__

    # Get resolution:edge lenght in m dict
    resolutions = {}
    for i in range(cell_min_res, cell_max_res, 1):
        resolutions[i] = h3.edge_length(resolution=i, unit='m')

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
    print(f"Start filling raster extent with h3 indexes at resolution {resolution}")
    df = pd.DataFrame({'cell_id': list(h3.polyfill_geojson(extent['features'][0]["geometry"], res=resolution))})

    # Get raster values for each cell_id
    print(f"Start getting raster values for hexes at resolution {resolution}")
    df[value_name] = df['cell_id'].apply(lambda x: raster_band_array[rs.index(h3.h3_to_geo(x)[1], h3.h3_to_geo(x)[0])])

    # Drop nodata
    df = df[df[value_name] != rs.nodata]

    return df


def vector_to_h3(vector_path, value_name, resolution, extent=None, layer=None):
    """Load raster values into h3 dggs cells

    Parameters:
    vector_path (string): path to vector file compatible with Geopandas for uploading
    value_name (string): vector attribute name to load into dggs cells
    resolution (integer): h3 cell resolution to use
    extent (list): extent as array of 2 lon lat pairs to get raster values for
    layer (string): vector layer name if geopackage is used

    Returns:
    Pandas dataframe
   """

    # Open vector to geodataframe
    gdf = gpd.read_file(vector_path, layer)

    # Get extent to fill with h3 hexes
    if extent:
        extent = gpd.GeoSeries(box(extent[0], extent[1], extent[2], extent[3])).__geo_interface__
    else:
        extent = gpd.GeoSeries(
            box(gdf['geometry'].total_bounds[0], gdf['geometry'].total_bounds[1], gdf['geometry'].total_bounds[2],
                gdf['geometry'].total_bounds[3])).__geo_interface__

    # Create dataframe with cell_ids from extent with given resolution
    print(f"Start filling raster extent with h3 indexes at resolution {resolution}")
    h3_gdf = gpd.GeoDataFrame({'cell_id': list(h3.polyfill_geojson(extent['features'][0]["geometry"], res=resolution))})

    # Get hex centroids for points
    h3_gdf['geometry'] = h3_gdf['cell_id'].apply(lambda x: Point(h3.h3_to_geo(x)[1], h3.h3_to_geo(x)[0]))
    hex_gdf = h3_gdf.set_crs('epsg:4326')

    # Spatial join hex centroids with gdf
    vector_h3 = gpd.sjoin(hex_gdf, gdf)

    # Drop unnecessary fields
    vector_h3 = vector_h3[['cell_id', value_name]]

    return vector_h3


def cell_h3_downsampling(df, cell_id_col, metric_col, coarse_resolution, metric_type):
    """Aggregates a given attribute in h3 cell to a given coarser resolution level

    Parameters:
    df (pandas dataframe): dataframe with s2 ids and attributes for aggregation
    cell_id_col (string): name of s2 id column
    metric_col (string): name of a column for aggreagation
    coarse_resolution (integer): Coarser s2 resoluiton for aggregation
    metric_type (string): attribute type (numerical, categorical)
    Returns:
    Pandas dataframe
   """

    df_coarse = df.copy()
    coarse_id_col = 'cell_id_{}'.format(coarse_resolution)
    df_coarse[coarse_id_col] = df_coarse[cell_id_col].apply(lambda x: h3.h3_to_parent(x, coarse_resolution))

    if metric_type == 'numeric':
        dfc = df_coarse.groupby(coarse_id_col)[[metric_col]].mean().reset_index()
    elif metric_type == 'categorical':
        dfc = df_coarse.groupby([coarse_id_col, metric_col]).agg(count=(metric_col, 'count')).reset_index().sort_values(
            by=[coarse_id_col, metric_col, 'count']).groupby(coarse_id_col, as_index=False, sort=False).first()
        dfc.drop('count', axis=1, inplace=True)
    dfc.columns = [cell_id_col, metric_col]
    return dfc



