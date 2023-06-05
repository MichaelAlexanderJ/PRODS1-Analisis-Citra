 # -*- coding: utf-8 -*-
"""
Created on Fri May 26 22:31:47 2023

@author: brian
"""

import os
from glob import glob
import geopandas as gpd
import fiona
import earthpy.spatial as es
import rioxarray as rxr
import matplotlib.pyplot as plt
import collections
import earthpy.plot as ep
import xarray as xr
from PIL import Image ,ImageOps, ImageEnhance
import numpy
import regex
import pandas
from bs4 import BeautifulSoup
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkt
import rasterio
from rasterio.enums import Resampling



def parseMetadata(filePath):
    with open(filePath) as file:
        soup = BeautifulSoup(file, features="lxml-xml")
        time=soup.find("PRODUCT_START_TIME")
        yearMonthString=time.string[:7].replace("-","")
        return yearMonthString

def createNormalizedPNG(outputPath,bandTuple,mode="RANGENORMALIZED"):
    resultImage = Image.new("I;16",tuple(reversed(bandTuple[0].shape)))
    resultArray=numpy.array(resultImage)
    resultArray=resultArray.astype(int)
    resultArray[:]=bandTuple[0].to_numpy()
    if mode=="RANGENORMALIZED":
        resultArray=resultArray/10000.0*numpy.iinfo("uint16").max
    elif mode=="MAXNORMALIZED":
        maxVal=resultArray.max()*1.0
        resultArray=resultArray/maxVal*numpy.iinfo("uint16").max 
    resultArray=resultArray.astype(numpy.uint16)
    resultImage=Image.fromarray(resultArray)
    resultImage.save(os.path.join(outputPath),"PNG",lossless=True)

def createRGB_R10(outputPath,R10Bands,mode="RANGENORMALIZED"):
    resultImage = Image.new("RGB",tuple(reversed(R10Bands[0][0].shape)))
    resultArray=numpy.array(resultImage)
    resultArray=resultArray.astype(int)
    resultArray[:,:,0]=R10Bands[2][0].to_numpy()
    resultArray[:,:,1]=R10Bands[1][0].to_numpy()
    resultArray[:,:,2]=R10Bands[0][0].to_numpy()
    if mode=="RANGENORMALIZED":
        resultArray=resultArray/10000.0*255
    elif mode=="MAXNORMALIZED":
        maxVal=resultArray.max()*1.0
        resultArray=resultArray/maxVal*255
    resultArray=resultArray.astype(numpy.uint8)
    resultImage=Image.fromarray(resultArray)
    resultImage.save(os.path.join(outputPath),"PNG",lossless=True)



def createMoistureIndex(outputPath,R20Bands):
    global rata
    global df_rata
    
    resultImageMI = Image.new("RGB",tuple(reversed(R20Bands[0][0].shape)))
    resultMI=numpy.array(resultImageMI)
    resultMI=resultMI.astype(int)
    
    band8a = R20Bands[9][0].to_numpy().astype(int)
    band11 = R20Bands[7][0].to_numpy().astype(int)
    
    mi1 = band8a-band11
    mi2 = band8a+band11
    
    mi = mi1/mi2
    rata = hitung_rata(mi)
    title = regex.search(r"\w+\s\w{4}\\(\w+\s\w+)",outputPath).group(1)
    
    data = {'Nama Tempat':title ,'Rata-Rata': rata}
    df_rata = df_rata.append(data,ignore_index=True)

    os.makedirs(outputPath,exist_ok=True)
    
    plot_MI(mi,outputPath,title)
    box_plot(mi,outputPath,title)
    
    


def box_plot(array,outplt,title):
    dtlist =[]
    for i in numpy.nditer(array):
        if not numpy.isnan(i):
            dtlist.append(i)
    plt.boxplot(dtlist)
    plt.title('Box Plot ' + title)
    plt.savefig(os.path.join(outplt,"boxplot "+title+ ".PNG"))



def plot_MI(mi,outputPath,title):
    plt.imshow(mi, cmap='Spectral', vmin=-1, vmax=1)
    plt.colorbar(label='Moisture Index')
    plt.title('Moisture Index ' + title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig(os.path.join(outputPath,"Mi "+title+".PNG"))
    plt.show()


def hitung_rata(array):
    count = 0
    nilai = 0
    for rows in array:
        for elemen in rows:
            if not numpy.isnan(elemen):
                nilai += elemen
                count += 1
    if count == 0:
        rata = float('nan')
    else:
        rata = nilai / count
    return rata


def create_allBandBoxPlot(dictAllband,outpath):
    allkeys=sorted(dictAllband.keys())
    allBand=[]
    for namaBand in allkeys:
        non_zero_data = dictAllband[namaBand].to_numpy()
        non_zero_data = non_zero_data[non_zero_data != 0]
        allBand.append(non_zero_data)

    
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.boxplot(allBand, labels=allkeys,sym='gx')
    ax.set_title('Sentinel 2')
    ax.title.set_fontsize(25)
    ax.set_xlabel('Bands')
    ax.xaxis.label.set_fontsize(20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(outpath,'boxplotBands.png'), dpi=300) 
    plt.show()
    
    
global df_rata
df_rata = pandas.DataFrame(columns=["Nama Tempat","Rata-Rata"])

outputFolder = os.path.join("E:\\croplanjutan","hasil crop")
geojsonPath = os.path.join("E:\\croplanjutan\\geojson","Sungaiii Karawang.geojson")
geometryLokasi = gpd.read_file(geojsonPath)



inputFolder=[
        os.path.join("E:\\croplanjutan\\band\\sungaiKarawang\\S2A_MSIL2A_20220718T025541_N0400_R032_T48MYU_20220718T081000.SAFE")
    ]

targetSubDirSet=set(["R10m", "R20m", "R60m"])

dictAllband = {}

for curInputFolder in inputFolder:
    granulFolders = []
    yearMonthString=""
    new_width,new_height=0,0
    for root, subdirs,files in os.walk(curInputFolder):
        subDirsAsSet = set(subdirs)
        
        if subDirsAsSet==targetSubDirSet:
            for curResFolder in subdirs:
                granulFolders.append(os.path.join(root,curResFolder))
            
        for curFile in files:
            if curFile=="MTD_MSIL2A.xml":
                yearMonthString=parseMetadata(os.path.join(root,curFile))
                
    for resFolder in granulFolders:
        jp2Files = glob(os.path.join(resFolder,"*_B??*.jp2"))
        crs = es.crs_check(jp2Files[0])
        targetGeometries = geometryLokasi.to_crs(crs)
        resolutionName = resFolder[-3:]
        
        curGeoOutputFolder = os.path.join(outputFolder,geojsonPath[geojsonPath.rindex(os.path.sep)+1:geojsonPath.rindex(".")],yearMonthString)
        curGeoOutputFolderCSV = os.path.join(curGeoOutputFolder,"CSV")
        curGeoOutputFolderRaster = os.path.join(curGeoOutputFolder,"Raster")
        curGeoOutputFolderRasterResize = os.path.join(curGeoOutputFolderRaster,"Resize")
        cropBands=[]
        cropBandsResize = []
        for i,bandFile in enumerate(jp2Files):
            bandName = regex.search(r"B\d+\w+_",bandFile)[0][:-1]
            notInBoundary=False
            with rxr.open_rasterio(bandFile, "r+") as dataset:
                dataset.rio.write_nodata(65535 ,inplace=True)
                try:
                    curBand=dataset.rio.clip(targetGeometries.geometry,from_disk=True).squeeze() 
                except Exception as e:
                    print(e)
                    notInBoundary=True 
                if notInBoundary:
                    break
                
                os.makedirs(curGeoOutputFolderCSV,exist_ok=True)
                os.makedirs(curGeoOutputFolderRaster,exist_ok=True)
                os.makedirs(curGeoOutputFolderRasterResize,exist_ok=True)
                curBand.rio.to_raster(os.path.join(curGeoOutputFolderRaster,f"{resolutionName}_{bandName}_RAWVALUE.TIFF"))
                
                if "R10" in resFolder:
                    if bandName not in dictAllband.keys():
                        dictAllband[bandName]=curBand
                        
                        
                # melakukan resize untuk radius 20m
                if "R20m" in resFolder:
                    resized_band = curBand.rio.reproject(curBand.rio.crs,shape=(new_width, new_height),resampling=Resampling.nearest)
                    resized_band.rio.to_raster(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_RAWVALUE.TIFF"))
                    if bandName not in dictAllband.keys():
                        dictAllband[bandName]=resized_band
                    
                # melakukan resize untuk radius 20m
                if "R60m" in resFolder:
                    resized_band = curBand.rio.reproject(curBand.rio.crs,shape=(new_width, new_height),resampling=Resampling.nearest)
                    resized_band.rio.to_raster(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_RAWVALUE.TIFF"))
                    if bandName not in dictAllband.keys():
                        dictAllband[bandName]=resized_band
                
            bandTuple=(curBand,bandName)
            createNormalizedPNG(os.path.join(curGeoOutputFolderRaster,f"{resolutionName}_{bandName}_RANGENORMALIZED.png"),bandTuple)
            createNormalizedPNG(os.path.join(curGeoOutputFolderRaster,f"{resolutionName}_{bandName}_MAXNORMALIZED.png"),bandTuple,mode="MAXNORMALIZED")
            cropBands.append(bandTuple)
           
            if "R20m" in resFolder:
                bandTupleResize=(resized_band,bandName)
                createNormalizedPNG(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_RANGENORMALIZED.png"),bandTupleResize)
                createNormalizedPNG(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_MAXNORMALIZED.png"),bandTupleResize,mode="MAXNORMALIZED")
                cropBandsResize.append(bandTupleResize)
                
            if "R60m" in resFolder:
                bandTupleResize=(resized_band,bandName)
                createNormalizedPNG(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_RANGENORMALIZED.png"),bandTupleResize)
                createNormalizedPNG(os.path.join(curGeoOutputFolderRasterResize,f"{resolutionName}_{bandName}_Resize_MAXNORMALIZED.png"),bandTupleResize,mode="MAXNORMALIZED")
                cropBandsResize.append(bandTupleResize)
            
                
        if notInBoundary:
            continue
        imageShape=cropBands[0][0].shape
        xIndex=[]
        yIndex=[]
        curX=0
        curY=0
        for i in range(0,imageShape[0]*imageShape[1]):
            xIndex.append(curX)
            yIndex.append(curY)
            curX=curX+1
            if curX==imageShape[1]:
                curX=0
                curY=curY+1
        resultDF=pandas.DataFrame()
        resultDF["x"]=pandas.Series(xIndex)
        resultDF["y"]=pandas.Series(yIndex)
        for bandData,bandName in cropBands:
            resultDF[bandName]= pandas.Series(bandData.to_numpy().flatten())
            
        resultDF.to_csv(os.path.join(curGeoOutputFolderCSV,f"{resolutionName}.csv"),index=False)
                
        
        if "R10m" in resFolder: ## bikin RGB image dari R10
            new_width,new_height=(cropBands[0][0].to_numpy().shape[0]),(cropBands[0][0].to_numpy().shape[1])
            createRGB_R10(os.path.join(curGeoOutputFolderRaster,"RGB10m.png"),cropBands)
            createRGB_R10(os.path.join(curGeoOutputFolderRaster,"RGB10m_MAXNORMALIZED.png"),cropBands,mode="MAXNORMALIZED")
        
        if "R20m" in resFolder: ## bikin MI image dari R20
            createMoistureIndex(os.path.join(os.path.join(curGeoOutputFolderRaster,"MI")),cropBands)
            
create_allBandBoxPlot(dictAllband,curGeoOutputFolderRaster)
