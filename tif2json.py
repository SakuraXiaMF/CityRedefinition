import json
import os
import gdal

tiffolder = "./TifFolder/"
chinatifpath = "osm_china_roadLength.tif"

def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res

filfiles = getsuffix(tiffolder,".tif")
feature_dicts = []

for tifpath in filfiles:
    print(tifpath)
    # tifpath = "./TifFolder/antarctica-latest-free_roadLength.tif"    # 1 代表 成功 0 代表失败 angola 1 benin 1 botswana 1  burkina 1 antarctica
    dataset = gdal.Open(tifpath)
    im_width = dataset.RasterXSize  # 栅格矩阵的列数
    im_height = dataset.RasterYSize  # 栅格矩阵的行数
    geo_transform = dataset.GetGeoTransform()
    origin_x = geo_transform[0]
    origin_y = geo_transform[3]
    pixel_width = geo_transform[1]
    pixel_height = geo_transform[5]
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 获取数据
    [rows, cols] = im_data.shape
    for i in range(rows):
        print("转换第 %d 行" % i)
        for j in range(cols):
            center_x = origin_x + j * pixel_width + pixel_width / 2
            center_y = origin_y + i * pixel_height + pixel_height / 2
            center_value = im_data[i, j]
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
result_dict = {
    "type": "FeatureCollection",
    "features": feature_dicts
}

with open("result.json", "w") as f:
    json.dump(result_dict, f)
    print("转换完成……")
