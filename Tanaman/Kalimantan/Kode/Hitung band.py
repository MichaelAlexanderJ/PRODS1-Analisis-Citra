# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:07:07 2023

@author: jerem
"""

import rasterio

# Membuka file .jp2
file_path = 'D:/Semester 8/PRODS1/Karawang/jp2 karawang cropped.jp2'
#file_path = 'D:/Semester 8/PRODS1/Minggu 8/transfer/Arjasari/Arjasari/202108/Raster/B01_R60m_MAXNORMALIZED.png'
dataset = rasterio.open(file_path)

# Mendapatkan jumlah band
num_bands = dataset.count

# Loop untuk menghitung nilai band 1 hingga 13
for band in range(1, num_bands + 1):
    # Membaca nilai piksel pada setiap band
    band_data = dataset.read(band)

    # Lakukan operasi atau pemrosesan yang diinginkan pada nilai piksel
    # Contoh: Menampilkan nilai piksel dari setiap band
    print("Band", band, "values:")
    print(band_data)

# Tutup dataset
dataset.close()

import rasterio

# Membuka file .jp2
file_path = 'D:/Semester 8/PRODS1/Karawang/jp2 karawang cropped.jp2'

try:
    # Open the file in read mode
    with rasterio.open(file_path) as dataset:
        # Mendapatkan jumlah band
        num_bands = dataset.count

        # List to store band data
        band_data_list = []

        # Loop untuk menghitung nilai band 1 hingga num_bands
        for band in range(1, num_bands + 1):
            # Membaca nilai piksel pada setiap band
            band_data = dataset.read(band)

            # Menambahkan band data ke dalam list
            band_data_list.append(band_data)

            # Contoh: Menampilkan nilai piksel dari setiap band
            print("Band", band, "values:")
            print(band_data)

        # Menampilkan semua band data
        print("All band data:")
        for i, band_data in enumerate(band_data_list, start=1):
            print("Band", i, "values:")
            print(band_data)

except rasterio.errors.RasterioIOError as e:
    # Handle any potential IO errors
    print("Error: ", e)

import rasterio

# Membuka file .jp2
file_path = 'D:/Semester 8/PRODS1/Karawang/jp2 karawang cropped.jp2'

try:
    # Open the file in read mode
    with rasterio.open(file_path) as dataset:
        # Mendapatkan jumlah band
        num_bands = dataset.count

        # Loop untuk menghitung nilai band 1 hingga num_bands
        for band in range(1, num_bands + 1):
            # Membaca nilai piksel pada setiap band
            band_data = dataset.read(band)

            # Menyimpan band data dalam variabel yang sesuai
            variable_name = f"band_data{band}"
            globals()[variable_name] = band_data

            # Contoh: Menampilkan nilai piksel dari setiap band
            print("Band", band, "values:")
            print(band_data)

        # Contoh: Mengakses band_data1 dan band_data2
        print("Accessing band_data1 and band_data2:")

except rasterio.errors.RasterioIOError as e:
    # Handle any potential IO errors
    print("Error: ", e)
