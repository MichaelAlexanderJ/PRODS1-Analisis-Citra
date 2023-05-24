# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:41:24 2023

@author: kris
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

global df_rata
df_rata = pandas.DataFrame(columns=["Nama Tempat","Rata-Rata"])

def resize(inputArray, newShape):
    inputImage=Image.new("I",tuple(reversed(inputArray.shape)))
    inputArray=numpy.array(inputImage)
    inputArray[:]=test1[:]
    inputImage=resultImage=Image.fromarray(inputArray,mode="I")
    newImage=inputImage.resize(tuple(reversed(newShape)),resample=Image.NEAREST)
    return numpy.array(newImage)

def parseMetadata(filePath):
    with open(filePath) as file:
        soup = BeautifulSoup(file, features="lxml-xml")
        time=soup.find("PRODUCT_START_TIME")
        yearMonthString=time.string[:7].replace("-","")
        return yearMonthString

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


def createMoistureIndex(outputPath,R20Bands,mode="RANGENORMALIZED"):
    # resultImageNatural = Image.new("RGB",tuple(reversed(R20Bands[0][0].shape)))
    # resultArray=numpy.array(resultImageNatural)
    # resultArray=resultArray.astype(int)
    
    # resultArray[:,:,0]=R20Bands[2][0].to_numpy()
    # resultArray[:,:,1]=R20Bands[1][0].to_numpy()
    # resultArray[:,:,2]=R20Bands[0][0].to_numpy()
    # maxVal=resultArray.max()*1.0
    # resultArray=resultArray/maxVal*255
    # resultArray=resultArray.astype(numpy.uint8)
    # resultImageNatural=Image.fromarray(resultArray)
    
    # # meningkatkat kontras 
    # enhancer = ImageEnhance.Contrast(resultImageNatural)
    # gamma_img = enhancer.enhance(2.0)
    
    
    # IMG Moisture Index
    global resultImageMI
    resultImageMI = Image.new("RGB",tuple(reversed(R20Bands[0][0].shape)))
    global resultMI
    resultMI=numpy.array(resultImageMI)
    resultMI=resultMI.astype(int)
    global mi
    global mi1
    global mi2
    mi1 = R20Bands[8][0].to_numpy()-R20Bands[6][0].to_numpy()
    mi2 = R20Bands[8][0].to_numpy()+R20Bands[6][0].to_numpy()
    global band8a
    band8a = R20Bands[8][0].to_numpy()
    global band11
    band11 = R20Bands[6][0].to_numpy()
    mi = mi1/mi2
    # mi = numpy.nan_to_num(mi,nan=0)
    global mindmi
    mindmi = (((mi - numpy.nanmin(mi))) / (numpy.nanmax(mi))) * 2 - 1
    colom = outputPath.split("\\")
    
    
    global rata
    rata = hitung_rata(mindmi)
    
    data = {'Nama Tempat':colom[5] , 'Rata-Rata': rata}
    # mindmi = (mindmi + 1)/2
    global df_rata
    df_rata = df_rata.append(data,ignore_index=True)
    
    # global G
    # global R
    
    # R = numpy.full(R20Bands[0][0].shape,0) - mindmi
    # R = (R - numpy.nanmin(R))/(0-numpy.nanmin(R))
    
    # G = numpy.full(R20Bands[0][0].shape,0) - mindmi
    # G = G + R
    # G = (G + 1)/2
    os.makedirs(outplt + '\\' + colom[5],exist_ok=True)
    plt.imshow(mindmi, cmap='Spectral', vmin=-1, vmax=1)
    plt.colorbar(label='Moisture Index')
    plt.title('Moisture Index Visualization')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig(os.path.join(outplt))
    plt.show()
    
    # resultMI[...,0]= R * 255
    # resultMI[...,1]= G * 255
    # resultMI[...,2]= mindmi * 255
    
    # # resultMI = numpy.interp(resultMI, (resultMI.min(), resultMI.max()), (0, 255))
    
    # resultMI=resultMI.astype(numpy.uint8)
    # resultImageMI=Image.fromarray(resultMI)
    
    # # Melakukan penggabungan gambar
    # global ress
    # ress=Image.blend(gamma_img,resultImageMI,alpha=0.5)
    # global resized_img
     
    # #  mengubah ukuran gambar
    # scale_factor = 2
    # resized_img = ress.resize((ress.width*scale_factor, ress.height*scale_factor), resample=Image.BICUBIC)
    
    # enhancer = ImageEnhance.Contrast(resized_img)
    # gamma_img = enhancer.enhance(2.0)
    # gamma_img.save(os.path.join(outputPath),"PNG",lossless=True)
    


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



outplt = os.path.join("E:\\","outplot")
outputFolder=os.path.join("E:\\","Sentinel Data ADM4")
admBoundariesPath = os.path.join("E:\\","gadm41_IDN.gpkg")
layers = fiona.listlayers(admBoundariesPath) # ngecek ada layer apa aja

admBoundaries=gpd.read_file(admBoundariesPath,layer="ADM_ADM_4") ## jgn select semua..bkln crash
jawaBaratADM4=admBoundaries[admBoundaries["NAME_1"]=="Jawa Barat"]
bandungADM4=jawaBaratADM4[jawaBaratADM4["NAME_2"].str.contains("bandung",case=False)]
kebonJerukADM4=jawaBaratADM4[jawaBaratADM4["NAME_4"].str.contains("Arjasari",case=False)]
andirADM4=jawaBaratADM4[jawaBaratADM4["NAME_3"]=="Andir"]
kotabogorAMD4=jawaBaratADM4[jawaBaratADM4["NAME_2"].str.contains("bogor",case=False)]
admBoundaries=gpd.read_file(admBoundariesPath,layer="ADM_ADM_2")
bandungADM2=admBoundaries[admBoundaries["NAME_2"]=="Bandung"]

# opGeometryDF=jawaBaratADM4[jawaBaratADM4["NAME_4"].str.contains("kebon jeruk",case=False)]
# opGeometryDF=opGeometryDF.to_crs(sentinel_crs)
# 
# test kebon jeruk aja
# targetGeometries=kebonJerukADM4

# #test bandung
targetGeometries=bandungADM4

# #test bogor
targetGeometries=kotabogorAMD4  

inputFolder=[
      # os.path.join("E:\\","S2B_MSIL2A_20210827T025539_N0301_R032_T48MYT_20210827T060508.SAFE"),
      os.path.join("E:\\","S2A_MSIL2A_20200827T025551_N0214_R032_T48MYT_20200827T064850.SAFE")
  ]



errorList=[]

targetSubDirSet=set(["R10m", "R20m", "R60m"] )
for curInputFolder in inputFolder:
    #cari folder granulenya & metadatanya
    
    granuleFolders=[]
    yearMonthString=""
    for root, subdirs, files in os.walk(curInputFolder):
        subDirAsSet=set(subdirs)
        if subDirAsSet==targetSubDirSet:
            for curResFolder in subdirs:
                granuleFolders.append(os.path.join(root,curResFolder))
   
        for curFile in files:
            if curFile=="MTD_MSIL2A.xml":
                yearMonthString=parseMetadata(os.path.join(root,curFile))
 

    for resFolder in granuleFolders:
        jp2Files=glob(os.path.join(resFolder,"*_B??*.jp2"))
        crs= es.crs_check(jp2Files[0])
        targetGeometries=targetGeometries.to_crs(crs)
        resolutionName=resFolder[resFolder.rindex(os.path.sep)+1:]

        for i,curGeomRow in targetGeometries.iterrows():
            curGeomOutputFolder=os.path.join(outputFolder,curGeomRow["NAME_1"],curGeomRow["NAME_2"],curGeomRow["NAME_3"],curGeomRow["NAME_4"],yearMonthString)
            curGeomOutputFolderCSV=os.path.join(curGeomOutputFolder,"CSV")
            curGeomOutputFolderRaster=os.path.join(curGeomOutputFolder,"Raster")
            os.makedirs(curGeomOutputFolder,exist_ok=True)

            croppedBands=[]
            for i,bandFile in enumerate(jp2Files):
                bandName=regex.search(r"_B.+?[_\.]+",bandFile)[0][1:-1]  
                notInBoundary=False
                with rxr.open_rasterio(bandFile, "r+") as dataset:
                    dataset.rio.write_nodata(65535 ,inplace=True)
                    try:
                        curBand=dataset.rio.clip(curGeomRow.geometry.geoms,from_disk=True).squeeze() 
                    except Exception as e:
                        print(e)
                        notInBoundary=True 
                          
                        errorRow=curGeomRow["NAME_1"]+":"+curGeomRow["NAME_2"]+":"+curGeomRow["NAME_3"]
                        errorList.append(errorRow)
                    if notInBoundary:
                        break
                    
                    #agak ngulang2 sih...
                    os.makedirs(curGeomOutputFolderCSV,exist_ok=True)
                    os.makedirs(curGeomOutputFolderRaster,exist_ok=True)
                    curBand.rio.to_raster(os.path.join(curGeomOutputFolderRaster,f"{bandName}_{resolutionName}_RAWVALUE.TIFF"))
                bandTuple=(curBand,bandName)
                createNormalizedPNG(os.path.join(curGeomOutputFolderRaster,f"{bandName}_{resolutionName}_RANGENORMALIZED.png"),bandTuple)
                createNormalizedPNG(os.path.join(curGeomOutputFolderRaster,f"{bandName}_{resolutionName}_MAXNORMALIZED.png"),bandTuple,mode="MAXNORMALIZED")
                croppedBands.append(bandTuple)
            if notInBoundary:
                continue
            imageShape=croppedBands[0][0].shape
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
            for bandData,bandName in croppedBands:
                resultDF[bandName]=pandas.Series(bandData.to_numpy().flatten())
            resultDF.to_csv(os.path.join(curGeomOutputFolderCSV,f"{resolutionName}.csv"),index=False)
            

    
            # if "R10m" in resFolder: ## bikin RGB image dari R10
            #     createRGB_R10(os.path.join(curGeomOutputFolderRaster,"RGB10m.png"),croppedBands)
            #     createRGB_R10(os.path.join(curGeomOutputFolderRaster,"RGB10m_MAXNORMALIZED.png"),croppedBands,mode="MAXNORMALIZED")
                          
            if "R20m" in resFolder: ## bikin MI image dari R20
                # R20Bands=croppedBands
                # break
                # print(os.path.join(curGeomOutputFolderRaster))
                createMoistureIndex(os.path.join(curGeomOutputFolderRaster,"PlotR20m.png"),croppedBands)
                # break
                # createMouMoistureIndex(os.path.join(curGeomOutputFolderRaster,"Moistureindex20m_MAXNORMALIZED.png"),croppedBands,mode="MAXNORMALIZED")

                          



# landsatPath=os.path.join("D:\landsat data","D:\landsat data\S2B_MSIL2A_20210827T025539_N0301_R032_T48MYT_20210827T060508\S2B_MSIL2A_20210827T025539_N0301_R032_T48MYT_20210827T060508.SAFE\GRANULE\L2A_T48MYT_A023367_20210827T031449\IMG_DATA\R10m")
# #landsatPath=os.path.join("D:\landsat data","D:\landsat data\S2B_MSIL1C_20180813T025539_N0206_R032_T48MYT_20180816T154649\S2B_MSIL1C_20180813T025539_N0206_R032_T48MYT_20180816T154649.SAFE\GRANULE\L1C_T48MYT_A007494_20180813T031633\IMG_DATA")
# jp2Files=glob(os.path.join(landsatPath,"*_B??*.jp2"))
# sentinel_crs = es.crs_check(jp2Files[0])
# #jp2Files.sort(key=lambda x:int(x[x.index("_B")+2:-4])) #gk perlu disort kynya, krn angka bandnya udah 2 digit




# admBoundaries=geopandas.read_file(admBoundariesPath,layer="ADM_ADM_4") ## jgn select semua..bkln crash
# jawaBaratADM4=admBoundaries[admBoundaries["NAME_1"]=="Jawa Barat"]

# admBoundaries=geopandas.read_file(admBoundariesPath,layer="ADM_ADM_2")
# bandung=admBoundaries[admBoundaries["NAME_2"]=="Bandung"]
# bandung.to_crs(sentinel_crs)




# # opGeometryDF=jawaBarat[jawaBarat["NAME_4"].str.contains("kebon jeruk",case=False)]
# # opGeometryDF=opGeometryDF.to_crs(sentinel_crs)

# opGeometryDF=jawaBaratADM4[jawaBaratADM4["NAME_4"].str.contains("kebon jeruk",case=False)]
# opGeometryDF=opGeometryDF.to_crs(sentinel_crs)
# currentGeometry=opGeometryDF.iloc[0].geometry
# #currentGeometry=bandung.iloc[0].geometry
# curBand=rxr.open_rasterio(jp2Files[0]).rio.clip(
#     currentGeometry,
#     from_disk=True).squeeze()




# croppedBands=[]
# for i,bandFile in enumerate(jp2Files):
#     with rxr.open_rasterio(bandFile, "r+") as dataset:
#         dataset.rio.write_nodata(-1 ,inplace=True)
#         curBand=dataset.rio.clip(currentGeometry,from_disk=True).squeeze() 
#     # curBand=rxr.open_rasterio(bandFile).rio.clip(
#     #     currentGeometry,
#     #     from_disk=True).squeeze() 
    
#     t=curBand.to_numpy()
    
#     bandName=regex.search(r"_B.+?[_\.]+",bandFile)[0][1:-1]
    
#     croppedBands.append((curBand,bandName))



# imageShape=croppedBands[0][0].shape
# flattenedArray=croppedBands[0][0].to_numpy().flatten()
# xIndex=[]
# yIndex=[]
# curX=0
# curY=0
# for i in range(0,imageShape[0]*imageShape[1]):
#     xIndex.append(curX)
#     yIndex.append(curY)
#     curX=curX+1
#     if curX==imageShape[1]:
#         curX=0
#         curY=curY+1

# resultDF=pandas.DataFrame()
# resultDF["x"]=pandas.Series(xIndex)
# resultDF["y"]=pandas.Series(yIndex)
# for bandData,bandName in croppedBands:
#     resultDF[bandName]=pandas.Series(bandData.to_numpy().flatten())








# counter = collections.Counter(curBand.to_numpy().flatten())


# croppedStack = xr.concat(croppedBands, dim="band")







# ep.plot_bands(croppedStack)

# ep.plot_rgb(croppedStack,
#             rgb=[3, 2, 1],
#             title="sentinel")




# plt.show()

