import os
from PIL import Image

# Path
# ganti path aja
TIFF_PATH = 'D:\\University_File\\Semester_8\\ProDS1\\project_analisis_citra\\PRODS1-Analisis-Citra\\Tanaman\\Kalimantan\\Resize'
OUTPUT_PATH = 'D:\\University_File\\Semester_8\\ProDS1\\project_analisis_citra\\PRODS1-Analisis-Citra\\Tanaman\\Kalimantan\\Hasil_resize'

# Size
NEW_WIDTH = 229
NEW_HEIGHT = 130
RESAMPLE = Image.NEAREST

# TIFF
TIFF = []
REGEX_PATTERN = r'*_B??*.TIFF'


def getTIFF():
    for tiff_name in os.listdir(TIFF_PATH):
        join = os.path.join(TIFF_PATH, tiff_name)
        tiff = Image.open(join)
        resize(tiff_name, tiff)


def resize(tiff_name, tiff_files):
    resize_image = tiff_files.resize(
        (NEW_WIDTH, NEW_HEIGHT),
        resample=RESAMPLE
    )
    saveTIFF(tiff_name, resize_image)


def saveTIFF(tiff_name, resize_image):
    save_nameEx = tiff_name + '.TIFF'
    save_path_name = os.path.join(OUTPUT_PATH, save_nameEx)
    resize_image.save(save_path_name)


def main():
    getTIFF()


if __name__ == '__main__':
    main()
