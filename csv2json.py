import json
import gdal
import numpy as np
import os



csvpath = "/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/China_ReDeFineCity.csv"
# dataset = gdal.Open("osm_china_roadLength.tif")
# im_width = dataset.RasterXSize  # 栅格矩阵的列数
# im_height = dataset.RasterYSize  # 栅格矩阵的行数
# geo_transform = dataset.GetGeoTransform()
# origin_x = geo_transform[0]
# origin_y = geo_transform[3]
# pixel_width = geo_transform[1]
# pixel_height = geo_transform[5]


def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res



# print(origin_x,origin_y,pixel_height,pixel_width)

# adfGeoTransform[0] /* top left x 左上角x坐标*/
# adfGeoTransform[1] /* w--e pixel resolution 东西方向上的像素分辨率*/
# adfGeoTransform[2] /* rotation, 0 if image is "north up" 如果北边朝上，地图的旋转角度*/
# adfGeoTransform[3] /* top left y 左上角y坐标*/
# adfGeoTransform[4] /* rotation, 0 if image is "north up" 如果北边朝上，地图的旋转角度*/
# adfGeoTransform[5] /* n-s pixel resolution 南北方向上的像素分辨率*/

origi = "RedefineCity.csv"#最开始的代码
files = getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/",".csv")
tiffiles = getsuffix("./TifFolder/",".tif")

files = sorted(files)
tiffiles = sorted(tiffiles)

input = []
feature_dicts = []
low_score = 100000

cnt=0
for file,tiffile in zip(files[cnt:],tiffiles[cnt:]):
    print(file,"\t",tiffile)
    tiffile = "/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/GRL_roadLength.tif"
    file = "/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/GRL_ReDeFineCity.csv"
    # tiffile = "osm_china_roadLength.tif"
    file = origi
    input = []
    dataset = gdal.Open(tiffile)
    im_width = dataset.RasterXSize  # 栅格矩阵的列数
    im_height = dataset.RasterYSize  # 栅格矩阵的行数
    geo_transform = dataset.GetGeoTransform()
    origin_x = geo_transform[0]
    origin_y = geo_transform[3]
    pixel_width = geo_transform[1]
    pixel_height = geo_transform[4]#geo_transform[5]
    print(origin_x,origin_y,pixel_height,pixel_width)
    fl  =0 
    for line in open(file, "r"):
        if fl ==0 :
            fl=1
            continue
        position = line.split('",')[0][1:]
        

        value = float(line.split('",')[1])

        if value < low_score and value !=0:
            low_score = value
        input.append([np.array(eval(position)), value])
    print( "读取完成……", len(input) )
    for line in input:
        center_x = 0
        center_y = 0
        # 网上的定位经纬度的代码
        # px = adfGeoTransform[0] + i * adfGeoTransform[1] + j * adfGeoTransform[2] 
        # py = adfGeoTransform[3] + i * adfGeoTransform[4] + j * adfGeoTransform[5]
        # 网上的代码到这---------
        for i in range(len(line[0])):
            x,y = line[0][i][0],line[0][i][1]# x 行 y 列
            px = origin_x + x*pixel_width+y*geo_transform[2]
            py = origin_y + x+pixel_height+y*geo_transform[5]
            center_x+=px
            center_y+=py
            # center_x += origin_x + y * pixel_width + pixel_width / 2
            # center_y += origin_y + x * pixel_height + pixel_height / 2
        center_x = center_x/len(line[0])
        center_y = center_y/len(line[0])
        center_value = (line[1] - low_score) * 100000
        # print(center_value,"Lat Long is ",center_y,center_x)
        feature_dict = {
            "type": "Feature",
            "properties": {
                "cum_conf": center_value,
                "latitude": center_y,
                "longitude": center_x
            },
            "geometry": {
                "type": "Point",
                "coordinates": [center_x, center_y]
            }
        }
        if center_value != 0:
            feature_dicts.append(feature_dict)
    break
result_dict = {
    "type": "FeatureCollection",
    "features": feature_dicts
}

with open("result.json", "w") as f:
    json.dump(result_dict, f)
    print("转换完成……")
