import json
import os

import gdal
import numpy as np



def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res

import json
import os

import gdal
import numpy as np

tiffile = "/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/GRL_roadLength.tif"



def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res
    




input = []
result_list = []
low_score = 100000

# filepath = "/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/GRL_ReDeFineCity.csv"#os.path.join("ReDeFineCity", "china.csv")
csvfiles = getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/",".csv")
TifFiles = getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/",".tif")

tiffilenames =[]
for tif in TifFiles:
    base = os.path.basename(tif)
    base = base.replace("__","_")
    name = base[:base.find("_roadLength.tif")]
    tiffilenames.append(name)

def check(name):
    # tmp = "/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/"+name+"_roadLength.tif"
    idx = -1
    # print(tmp)
    if name in tiffilenames:
        idx = tiffilenames.index(name)
    return idx
#  _ReDeFineCity    _roadLength
if __name__ == "__main__":
    for csvfile in csvfiles:
        base = os.path.basename(csvfile)
        base = base.replace("__","_")
        name = base[:base.find("_ReDeFineCity.csv")]
        if name in tiffilenames:
            idx = tiffilenames.index(name)
        else:
            print(name,"Fault,No TIF file")
            continue
        tiffile = TifFiles[idx]
        print(name,"True",tiffile)
        dataset =gdal.Open(tiffile) #gdal.Open("osm_china_roadLength.tif")
        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数
        geo_transform = dataset.GetGeoTransform()
        origin_x = geo_transform[0]
        origin_y = geo_transform[3]
        pixel_width = geo_transform[1]
        pixel_height = geo_transform[5]

        for line in open(csvfile, "r"):
            position = line.split('",')[0][1:]
            if 'SumLength' in position:
                continue
            value = float(line.split('",')[1])
            if value < low_score:
                low_score = value
            input.append([np.array(eval(position)), value])

        for line in input:
            center_x = 0
            center_y = 0
            for i in range(len(line[0])):
                x, y = line[0][i][0], line[0][i][1]
                center_x += origin_x + y * pixel_width + pixel_width / 2
                center_y += origin_y + x * pixel_height + pixel_height / 2
            center_x = center_x / len(line[0])
            center_y = center_y / len(line[0])
            center_value = (line[1] - low_score) * 100000
            if center_value != 0:
                result_list.append(str(center_value) + "," + str(center_y) + "," + str(center_x) + "\n")

    with open("result1.csv", "w") as f:
        f.write("cum_conf,latitude,longitude\n")
        f.writelines(result_list)
        result_list = []
    print("转换完成……")
    # print("手动写入")
    # path = "/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/ByAuto/"
    # autofiles = getsuffix(path,".csv")

    # import csv

    # for file in autofiles:
    #     print(file)
    #     with open(file,'r') as csvfile:
    #         reader = csv.reader(csvfile)
    #         rows = [str(row)[1:-1] for row in reader]
    #     with open("result1.csv","a") as f:
    #         f.writelines(rows)
    #         rows=[]
    # print("写入完成")