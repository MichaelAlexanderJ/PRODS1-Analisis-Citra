# -*- coding: utf-8 -*-
"""
Created on Sun May 21 13:46:27 2023

@author: jerem
"""

import rasterio
from rasterio.mask import mask
import geopandas as gpd
from pyproj import CRS


#           Gagal karena  ValueError: Input shapes do not overlap raster.
# # Baca file .jp2
# jp2_file = 'D:/Semester 8/PRODS1/Kalimantan/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158.SAFE/GRANULE/L2A_T50MMD_A019062_20190215T024740/IMG_DATA/R10m/T50MMD_20190215T022801_WVP_10m.jp2'
# dataset = rasterio.open(jp2_file)

# # Baca file .geojson
# geojson_file = 'D:/Semester 8/PRODS1/Kalimantan/geojson hutan sawit kalimantan.geojson'
# geojson = gpd.read_file(geojson_file)

# # Ambil geometri dari file .geojson
# geometry = geojson.geometry.values[0]

# # Ubah geometri menjadi format yang diterima oleh rasterio
# crop_shape = [geometry.__geo_interface__]

# # Lakukan transformasi CRS dari file .geojson ke CRS file .jp2
# target_crs = CRS.from_string(dataset.crs.to_string())
# geojson_crs = CRS.from_user_input(geojson.crs)
# transformed_geojson = geojson.to_crs(target_crs)

# # Lakukan cropping pada file .jp2
# cropped_image, cropped_transform = mask(dataset, crop_shape, crop=True)

# # Dapatkan metadata dari file .jp2 yang dicrop
# cropped_meta = dataset.meta.copy()
# cropped_meta.update({
#     'height': cropped_image.shape[1],
#     'width': cropped_image.shape[2],
#     'transform': cropped_transform,
#     'crs': target_crs
# })

# # Simpan file .jp2 yang dicrop dengan CRS yang sama
# output_file = 'D:/Semester 8/PRODS1/Kalimantan/crop1.jp2'
# with rasterio.open(output_file, 'w', **cropped_meta) as dest:
#     dest.write(cropped_image)

# # Update CRS file .geojson
# transformed_geojson.to_file('geojson hutan sawit kalimantan cropped.geojson', driver='GeoJSON')


#                       versi CRS tidak sama
import rasterio
import geopandas as gpd

# Open the raster image and get its CRS
jp2_file = 'D:/Semester 8/PRODS1/Kalimantan/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158.SAFE/GRANULE/L2A_T50MMD_A019062_20190215T024740/IMG_DATA/R10m/T50MMD_20190215T022801_WVP_10m.jp2'
dataset = rasterio.open(jp2_file)
raster_crs = dataset.crs

# Read the GeoJSON file and get its CRS
geojson_file = 'D:/Semester 8/PRODS1/Kalimantan/geojson hutan sawit kalimantan.geojson'
geojson = gpd.read_file(geojson_file)
geojson_crs = geojson.crs

# Print the CRS information
print("Raster CRS:", raster_crs)
print("GeoJSON CRS:", geojson_crs)


#                   Samakan versi CRS
import rasterio
import geopandas as gpd
from pyproj import CRS
from rasterio.warp import calculate_default_transform, reproject

# Baca file .jp2
jp2_file = 'D:/Semester 8/PRODS1/Kalimantan/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158.SAFE/GRANULE/L2A_T50MMD_A019062_20190215T024740/IMG_DATA/R10m/T50MMD_20190215T022801_WVP_10m.jp2'
dataset = rasterio.open(jp2_file)

# Baca file .geojson
geojson_file = 'D:/Semester 8/PRODS1/Kalimantan/geojson hutan sawit kalimantan.geojson'
geojson = gpd.read_file(geojson_file)

# Cek CRS
raster_crs = dataset.crs
geojson_crs = geojson.crs

if raster_crs != geojson_crs:
    # Reproyeksi GeoJSON ke CRS raster
    target_crs = raster_crs
    transformed_geojson = geojson.to_crs(target_crs)

    # Simpan GeoJSON yang telah direproyeksi
    transformed_geojson.to_file('D:/Semester 8/PRODS1/Kalimantan/geojson kalimantan transformed.geojson', driver='GeoJSON')

    # Baca ulang GeoJSON yang telah direproyeksi
    reprojected_geojson = gpd.read_file('D:/Semester 8/PRODS1/Kalimantan/geojson kalimantan transformed.geojson')

    # Reproyeksi raster ke CRS GeoJSON
    transform, width, height = calculate_default_transform(dataset.crs, target_crs, dataset.width, dataset.height, *dataset.bounds)
    reprojected_array = rasterio.features.geometry_mask(reprojected_geojson.geometry, out_shape=(height, width), transform=transform, invert=True)

    # Buat metadata baru untuk raster yang direproyeksi
    reprojected_meta = dataset.meta.copy()
    reprojected_meta.update({
        'crs': target_crs,
        'transform': transform,
        'width': width,
        'height': height
    })

    # Simpan raster yang direproyeksi
    with rasterio.open('D:/Semester 8/PRODS1/Kalimantan/jp2 kalimantan transformed.jp2', 'w', **reprojected_meta) as dst:
        dst.write(reprojected_array, 1)
else:
    print("CRS raster dan GeoJSON sudah sama.")


#           Cek CRS, CRS sudah sama
import rasterio
import geopandas as gpd

# Open the raster image and get its CRS
jp2_file = 'D:/Semester 8/PRODS1/Kalimantan/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158/S2A_MSIL2A_20190215T022801_N0211_R046_T50MMD_20190215T062158.SAFE/GRANULE/L2A_T50MMD_A019062_20190215T024740/IMG_DATA/R10m/T50MMD_20190215T022801_WVP_10m.jp2'
dataset = rasterio.open(jp2_file)
raster_crs = dataset.crs

# Read the GeoJSON file and get its CRS
geojson_file = 'D:/Semester 8/PRODS1/Kalimantan/geojson kalimantan transformed.geojson'
geojson = gpd.read_file(geojson_file)
geojson_crs = geojson.crs

# Print the CRS information
print("Raster CRS:", raster_crs)
print("GeoJSON CRS:", geojson_crs)



#           Cropping .jp2 berdasarkan .geojson
import rasterio
from rasterio.mask import mask
import geopandas as gpd

# Baca file .jp2
jp2_file = 'D:/Semester 8/PRODS1/Kalimantan/jp2 kalimantan transformed.jp2'
dataset = rasterio.open(jp2_file)

# Baca file .geojson
geojson_file = 'D:/Semester 8/PRODS1/Kalimantan/geojson kalimantan transformed.geojson'
geojson = gpd.read_file(geojson_file)

# Ambil geometri dari file .geojson
geometry = geojson.geometry.values[0]

# Ubah geometri menjadi format yang diterima oleh rasterio
crop_shape = [geometry.__geo_interface__]

# Lakukan cropping pada file .jp2
cropped_image, cropped_transform = mask(dataset, crop_shape, crop=True)

# Dapatkan metadata dari file .jp2 yang dicrop
cropped_meta = dataset.meta.copy()
cropped_meta.update({
    'height': cropped_image.shape[1],
    'width': cropped_image.shape[2],
    'transform': cropped_transform
})

# Simpan file .jp2 yang dicrop
output_file = 'D:/Semester 8/PRODS1/Kalimantan/jp2 kalimantan cropped.jp2'
with rasterio.open(output_file, 'w', **cropped_meta) as dest:
    dest.write(cropped_image)




