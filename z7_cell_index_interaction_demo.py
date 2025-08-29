# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-2025 - Alexander Kmoch
# Licenced under GNU AFFERO GENERAL PUBLIC LICENSE. Please consult the LICENCE
# file for details.
#
# Authors:
# Alexander Kmoch (alexander.kmoch@ut.ee)
# Wai Tik Chan (wai.tik.chan@ut.ee)

import geopandas as gpd
import numpy as np
import pandas as pd


def decode_z7_index(hex_string):
    """
    Decode a Z7 hexadecimal index (provided as hex string).
    
    Format:
    - First 4 bits: base cell number (0-11)
    - Remaining 60 bits: 20 groups of 3 bits each for resolution digits (0-6, 7 for beyond resolution)
    
    Args:
        hex_string (str): Hexadecimal string representing the Z7 cell index
        
    Returns:
        tuple: (base_cell, resolution_digits)
            - base_cell: integer 0-11
            - resolution_digits: list of integers (0-6 or 7), with 7 indicating beyond resolution
    """
    # Convert hex to binary string, ensuring we have full 64 bits
    binary = bin(int(hex_string, 16))[2:].zfill(64)
    
    # Extract base cell (first 4 bits)
    base_cell = int(binary[:4], 2)
    
    # Extract resolution digits (20 groups of 3 bits each)
    resolution_digits = []
    for i in range(20):  # 60 remaining bits = 20 groups of 3 bits
        start = 4 + (i * 3)  # Start after the first 4 bits
        value = int(binary[start:start + 3], 2)
        resolution_digits.append(value)
    
    return base_cell, resolution_digits


def z7_to_z7string(hex_string):
    """
    Get the Z7 string representation of a Z7 hexadecimal representation.

    Args:
        hex_string (str):  Z7 hexadecimal string represention of the cell index

    """
    base_cell, resolution_digits = decode_z7_index(hex_string)
    str_rep = [str(base_cell).zfill(2)]
    for digit in resolution_digits:
        status = False if digit == 7 else True
        if status:
            str_rep.append(str(digit))
    return "".join(str_rep)


def hex_to_int( hex_str):
    # From hex string to integer
    value = int(hex_str, 16)
    return value


def int_to_hex(value):
    # From integer back to hex string, maintaining 16 characters (64 bits)
    hex_back = f"{value:016x}"
    return hex_back


def get_z7_resolution(hex_str):
    """
    Get the resolution of a Z7 cell.

    Args:
        hex_str (str): Z7 hexadecimal string represention of the cell index
    """
    base_cell, resolution_digits = decode_z7_index(hex_str)    
    return len( list( filter(lambda x: x >= 0 and x < 7, resolution_digits) ) )


def get_z7string_resolution(z7str):
    """
    Get the resolution of a Z7 cell from its Z7_STRING representation.

    Args:
        z7str (str): Z7_STRING representation of the cell index
    """
    return len(z7str) - 2


def get_z7_local_pos(hex_str):
    """
    Get the local position of a cell within its parent cell.

    Args:
        hex_str (str): Z7 hexadecimal string represention of the cell index
    """
    z7_string = z7_to_z7string(hex_str)
    parent = z7_string[0:-1]
    local_pos = z7_string[-1:]
    is_center = True if local_pos == "0" else False
    return (parent, local_pos, is_center)


def get_z7string_local_pos(z7_string):
    """
    Get the local position of a cell within its parent cell.

    Args:
        hex_str (str): Z7 string represention of the cell index
    """
    parent = z7_string[0:-1]
    local_pos = z7_string[-1:]
    is_center = True if local_pos == "0" else False
    return (parent, local_pos, is_center)


def get_neighbours_by_index(gdf_idx, gpd_sindex):
    """
    Get the dataframe indices of the neighbouring cells of a given cell index.

    Args:
        gdf_idx (int): Index of a GeoDataframe row of the cell for which to find neighbours
        gpd_sindex (numpy.ndarray): Spatial index of the GeoDataFrame, pre-compute via
            `gpd_sindex = gdf.sindex.query(gdf.geometry, predicate="intersects")`
    """
    # Get positions where first row equals 0
    positions_1 = np.where(gpd_sindex[0] == gdf_idx)[0]
    # Get values from second row at those positions
    values_1 = gpd_sindex[1, positions_1]
    # and inverse to be sure
    positions_2 = np.where(gpd_sindex[1] == gdf_idx)[0]
    values_2 = gpd_sindex[1, positions_2]
    return set(values_1.tolist()).union(set(values_2.tolist()))


def get_neighbours_by_z7(z7_idx, gdf, gpd_sindex, z7_col='name'):
    """
    Get the Z7/Z7_STRING indices of the neighbouring cells of a given cell index.

    Args:
        z7_idx (str): Z7/Z7_STRING index of the cell for which to find neighbours
        gdf (GeoDataFrame): GeoDataFrame containing the cells
        gpd_sindex (numpy.ndarray): Spatial index of the GeoDataFrame, pre-compute via
            `gpd_sindex = gdf.sindex.query(gdf.geometry, predicate="intersects")`
        z7_col (str): Name of the column containing the Z7/Z7_STRING indices
    """
    idx_arr = gdf.loc[gdf[z7_col] == z7_idx].index.tolist()
    if len(idx_arr) <= 0:
        raise ValueError(f"{z7_idx} not in {z7_col}")
    idx = idx_arr[0]
    values_set = get_neighbours_by_index(gdf_idx=idx, gpd_sindex=gpd_sindex)
    values_set = values_set - set([idx])
    neighbour_indices = gdf.iloc[list(values_set)][z7_col].tolist()
    return neighbour_indices


# convenience function to apply to a GeoDataFrame
def apply_convert_z7(hex_str):
    z7_string = z7_to_z7string(hex_str)
    z7_res = get_z7_resolution(hex_str)
    parent, local_pos, is_center = get_z7_local_pos(hex_str)
    return pd.Series([z7_string, z7_res, parent, local_pos, is_center])
    

if __name__ == "__main__":
    grid_fname = '/Users/akmoch/Nextcloud/shared_work/ut_work/supervision/aleksandra_rammul_dggs_scale_metrics/data/alutag_igeo7_res_9.gpkg'
    gdf = gpd.read_file(grid_fname)

    gdf[['z7_string', 'z7_res', 'parent', 'local_pos', 'is_center' ]] = gdf['name'].apply(apply_convert_z7)

    # do in notebook only
    # gdf.explore()

    gpd_sindex = gdf.sindex.query(gdf.geometry, predicate="intersects")

    z7_hex_str = '0042aad3ffffffff'

    # lower right corner
    n = get_neighbours_by_z7(z7_idx=z7_hex_str,
                             gdf=gdf, gpd_sindex=gpd_sindex,
                             z7_col='name')
    
    z7_string = z7_to_z7string(z7_hex_str)
    print(f"""Z7 string representation of {z7_hex_str}: {z7_string}""")
    print(f"""Resolution of {z7_hex_str}: {get_z7_resolution(z7_hex_str)}""")
    
    print(f"""Neighbours of {z7_hex_str}: {n}""")

    parent, local_pos, is_center = get_z7_local_pos(z7_hex_str)
    print(f"""Parent (in z7_string format) cell of {z7_hex_str}: {parent}""")