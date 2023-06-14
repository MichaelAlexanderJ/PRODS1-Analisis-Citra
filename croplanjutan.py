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

    # ax.boxplot(allBand, labels=allkeys, sym='gx', medianprops=dict(color=median_color), boxprops=boxprops, whiskerprops=boxprops, capprops=boxprops)


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


def createPerBandBoxPlot(bands,outpath):
    allkeyss=sorted(bands.keys())
    allBand=[]
    for namaBand in allkeyss:
        non_zero_data = bands[namaBand]
        non_zero_data = non_zero_data[non_zero_data != 0]
        allBand.append(non_zero_data)

    # ax.boxplot(allBand, labels=allkeys, sym='gx', medianprops=dict(color=median_color), boxprops=boxprops, whiskerprops=boxprops, capprops=boxprops)

    fig, ax = plt.subplots(figsize=(20, 10))
    ax.boxplot(allBand, labels=allkeyss,sym='gx')
    ax.set_title('Sentinel 2')
    ax.title.set_fontsize(25)
    ax.set_xlabel('Bands')
    ax.xaxis.label.set_fontsize(20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(outpath,'boxplotBandPilihan.png'), dpi=300) 
    plt.show()
    

def replace_values(arr,minimal,masimal):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i, j] < minimal or arr[i, j] > masimal:
                arr[i, j] = 0
    return arr


def flood_fill(arr,boolean,x ,y, new):
    global count_panjang
    if x < 0 or x >= len(arr[0]) or y < 0 or y >= len(arr[0]):
        return

    if arr[x][y] == 0 or boolean[x][y] == True:
        return
    
    arr[x][y] = new
    boolean[x][y]= True
    count_panjang +=1
    flood_fill(arr,boolean,x+1, y, new)
    flood_fill(arr,boolean,x-1, y, new)
    flood_fill(arr,boolean,x, y+1, new)
    flood_fill(arr,boolean,x, y-1, new)
    


def flood(arrayBand):
    global index_Sungai,panjang_sungai,count_panjang,nilai_sungai
    global field,boolean
    
    field = numpy.copy(arrayBand)
    boolean = numpy.zeros_like(field,dtype=bool)
    panjang_sungai,count_panjang,nilai_sungai = 0,0,0
    index_Sungai = 1
    
    for i in range(len(field)):
        for j in range(len(field[0])):
            count_panjang = 0
            flood_fill(field, boolean,i ,j,index_Sungai)
            
            if (count_panjang > panjang_sungai):
                panjang_sungai = count_panjang 
                nilai_sungai = index_Sungai
            index_Sungai +=1
    
    
    field = numpy.array(field)
    field[field != nilai_sungai] = 0
    global indices
    indices = numpy.where(field == nilai_sungai)
    new_flid = numpy.zeros_like(arrayBand)
    for i in range(len(indices[0])):
        new_flid[indices[0][i]][indices[1][i]] = array1[indices[0][i]][indices[1][i]]

    return new_flid



global df_rata
df_rata = pandas.DataFrame(columns=["Nama Tempat","Rata-Rata"])

outputFolder = os.path.join("E:\\croplanjutan","hasil crop")
geojsonPath = os.path.join("E:\\croplanjutan\\geojson","Sungai Karawang.geojson")
geometryLokasi = gpd.read_file(geojsonPath)



inputFolder=[
        os.path.join("E:\\croplanjutan\\band\\sungaiKarawang\\S2A_MSIL2A_20220718T025541_N0400_R032_T48MYU_20220718T081000.SAFE")
    ]

targetSubDirSet=set(["R10m", "R20m", "R60m"])

dictAllband = {}

for curInputFolder in inputFolder:
    granulFolders = []
    yearMonthString=""
    new_width,new_height= 0,0
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
BandPilihan={}


    
def createArrayNormalizedPNG(array,outputPath,name,mode="RANGENORMALIZED"):
    global resultImage
    print(array.shape)
    outputPath = os.path.join(outputPath,"flood fill")
    os.makedirs(outputPath,exist_ok=True)
    resultImage = Image.new("I;16",tuple(reversed(array.shape)))
    resultArray=numpy.array(resultImage)
    resultArray=resultArray.astype(int)
    resultArray[:]=array
    if mode=="RANGENORMALIZED":
        resultArray=resultArray/10000.0*numpy.iinfo("uint16").max
    elif mode=="MAXNORMALIZED":
        maxVal=resultArray.max()*1.0
        resultArray=resultArray/maxVal*numpy.iinfo("uint16").max 
    resultArray=resultArray.astype(numpy.uint16)
    resultImage=Image.fromarray(resultArray)
    resultImage.save(os.path.join(outputPath, f"{name}.png"), "PNG", lossless=True)



array1 = numpy.copy(dictAllband['B8A'].to_numpy())
name='B8A'
minimalb8A = 1000
maksimalb8A = 3000

array1 = replace_values(array1,minimalb8A,maksimalb8A)

array1= flood(array1)

BandPilihan[name]=array1
# createPerBandBoxPlot(new)

createArrayNormalizedPNG(array1,curGeoOutputFolderRaster,name)

# createArrayNormalizedPNG(new,mode="MAXNORMALIZED")


array2 = numpy.copy(dictAllband['B08'].to_numpy())
name='B08'
minimalb08 = 1000
maksimalb08 = 3000

array2 = replace_values(array2,minimalb08,maksimalb08)

array2= flood(array2)
tuplePilihan = (array2,name)

BandPilihan[name]=array2

createArrayNormalizedPNG(array2,curGeoOutputFolderRaster,name)

array3 = numpy.copy(dictAllband['B09'].to_numpy())
name='B09'
minimalB09 = 2000
maksimalB09 = 3300

array3 = replace_values(array3,minimalB09,maksimalB09)

array3= flood(array3)
tuplePilihan = (array3,name)

BandPilihan[name]=array3

createArrayNormalizedPNG(array3,curGeoOutputFolderRaster,name)

array4 = numpy.copy(dictAllband['B12'].to_numpy())
name = 'B12'
minimalB12= 1000
maksimalB12 = 1800

array4 = replace_values(array4,minimalB12,maksimalB12)

array4= flood(array4)
tuplePilihan = (array4,name)

BandPilihan[name]=array4

createArrayNormalizedPNG(array4,curGeoOutputFolderRaster,name)

createPerBandBoxPlot(BandPilihan,curGeoOutputFolderRaster)















