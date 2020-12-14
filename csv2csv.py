import json
import os

import gdal
import numpy as np

dataset = gdal.Open("osm_china_roadLength.tif")
im_width = dataset.RasterXSize  # 栅格矩阵的列数
im_height = dataset.RasterYSize  # 栅格矩阵的行数
geo_transform = dataset.GetGeoTransform()
origin_x = geo_transform[0]
origin_y = geo_transform[3]
pixel_width = geo_transform[1]
pixel_height = geo_transform[5]

input = []
result_list = []
low_score = 100000
filepath = os.path.join("ReDeFineCity", "GRL_ReDeFineCity.csv")
for line in open(filepath, "r"):
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
        result_list.append(str(center_x) + "," + str(center_y) + "," + str(center_value) + "\n")

with open("result.csv", "w") as f:
    f.write("lng,lat,count\n")
    f.writelines(result_list)
    result_list = []
print("转换完成……")
