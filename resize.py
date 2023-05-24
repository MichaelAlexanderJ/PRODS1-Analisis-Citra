#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 13:20:06 2023

@author: vincentiusdaryl
"""

import cv2
import glymur
import numpy as np

def read_jp2_image(image_path):
    jp2 = glymur.Jp2k(image_path)
    return np.array(jp2[:])

def write_jp2_image(image, output_path):
    glymur.Jp2k(output_path, data=image)

def resize_images(image_paths, target_size):
    resized_images = []
    for image_path in image_paths:
        image = read_jp2_image(image_path)
        resized_image = cv2.resize(image, target_size, interpolation=cv2.INTER_AREA)
        resized_images.append(resized_image)
    return resized_images

# Ganti dengan path gambar Anda
image_paths = ['kalimantan_r10_crop.jp2', 'kalimantan_r20_crop.jp2', 'kalimantan_r20_crop.jp2']

# Ganti dengan ukuran pixel yang diinginkan (lebar, tinggi)
target_size = (1820, 1364)

resized_images = resize_images(image_paths, target_size)

# Menyimpan gambar yang telah diubah ukurannya
for i, resized_image in enumerate(resized_images):
    write_jp2_image(resized_image, f'gambar{i+1}_resized.jp2')
