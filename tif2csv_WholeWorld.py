#  -*-  coding:utf-8  -*-
import requests
from tqdm import tqdm
import json
import logging
import re
import numpy as np
from geopy.geocoders import Nominatim
import ssl
import ast
import sys
import csv
import time as tm
import sys
import datetime
import matplotlib.pyplot as plt
import math
import os
import re
import pandas as pd
from multiprocessing import Process, Pool, current_process
import random
from shapely.geometry import Point, Polygon
import shapely.wkt
from geopy.distance import geodesic
import MySQLdb
from sqlalchemy import create_engine
from IPython.display import clear_output
from datetime import datetime
import geopandas as gpd

import geopandas
import matplotlib.pyplot as plt
pd.set_option('max_colwidth',200)
# chinagpd = geopandas.read_file("./gis_osm_roads_free_1.shp")
# chinagpd.sample(2)

# 全球-重新定义城市

import os
from osgeo import gdal
from osgeo import ogr
import numpy as np
import matplotlib.pyplot as plt
import datetime






def getroute(order,sumlen):
    tmp = []
    global route
    for o in order:
        tmp.append(o)
    route.append([tmp,sumlen])
    
def Walking(x,y,ds_array,xlimit,ylimit,p,sumlen,limitout):
    global vis,pq,order;
    pq.append((x,y))
    tmpcnt=0
    fl = 0
    while pq:
        v = pq.pop(0)
        x,y=v[0],v[1]
        order.append([x,y])
        for i in range(len(stepx)):
            nx,ny = x+stepx[i],y+stepy[i]
            if nx>=xlimit or ny>=ylimit or nx<0 or ny < 0:
                continue
            tmp = ds_array[nx][ny]
            if vis[nx][ny]==0 and tmp - limitout > 1e-8 and tmp != 0 :#and sumlen+tmp<=p:
                sumlen+=tmp
                vis[nx][ny]=1
                if  p-sumlen <= 1e-8:
                    fl =1 
                    getroute(order,sumlen)
                    return 0;
                pq.append( (nx,ny) )
    if len(order)>0 and p-sumlen <=1e-5:
        getroute(order,sumlen)
    return 0;

    


def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res


pq = [ ]
order = [ ] 
sumlen=0
route = [ ]
stepx=[-1,  0,  1, -1, 1, -1, 0, 1]
stepy=[-1, -1, -1,  0, 0,  1, 1, 1]

# stepx=[  0, -1, 1, 0 ]
# stepy=[ -1,  0, 0, 1 ]
# vis = [ [0 for i in range(limit)] for i in range(limit) ]
    
    
if __name__ == "__main__":
    starttime = datetime.datetime.now()
    TifFiles = getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/",".tif")
    cnt=0
    # cnt = TifFiles.index("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/china_roadLength.tif")
    for filepath in TifFiles[cnt:]:
        # filepath = "/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/GRL_roadLength.tif"
        pq = []
        order = []
        sumlen=0
        route = [ ]
        print(filepath)
        base = os.path.basename(filepath)
        name = base[:base.find("_roadLength")]

        ds = gdal.Open(filepath)
        ds_array = ds.ReadAsArray()
    
        x = ds.RasterYSize  # 栅格矩阵的行数
        y = ds.RasterXSize  # 栅格矩阵的列数
        # limit =  ds.GetRasterBand(1).GetNoDataValue()
        limit = -float("inf")
        # if name == "JPN" or name == "CAN":
        print("limit is ",limit)
        # for i in range(0,x):
        #     # print(ds_array[i][0])
        #     for j in range(0,y):
        #         if ds_array[i][j] == limit :
        #             # print(ds_array[i][j])
        #             ds_array[i][j] = 0
        #             # print(ds_array[i][j])

        print("Square is ({} * {})".format(x,y))
        vis = [ [0 for i in range(y)] for i in range(x) ]
        data = []
        sumrow = 0
        excess = 0
        cnt  = 0 
        maxs,mins =-1,1000000000000000
        for rindx,row in enumerate(ds_array):
            # if name == "JPN" or name == "CAN":
                # sumrows = round(sum(row),10)
                # print(name,sumrows)
            # if sumrow==0:
            #     continue
            # sumrow +=round(sumrows,10)
            # maxs,mins = max(row),min(row)
            # if name =="JPN" or name == "CAN":
            #     excess = 0.00005
            # else:
            for cindx,cow in enumerate(row):
                if cow ==0 or cow == limit or cow == None :
                    continue
                maxs = max(maxs,cow)
                mins = min(mins,cow)
                cow = round(cow,12)
                sumrow+=cow
                # if name == "JPN" or name == "CAN":
                #     print(cow)
                data.append( [rindx,cindx,cow] )
        percent = 0#round(len(data)/(x*y),5)
        excess = max(excess,maxs-mins)# 阈值计算方法
        # p = ( sumrow/len(data) ) * min( 7 , len(data) )
        p = excess*20 #excess*(6)
        if sumrow<p:
            p=excess
        if p == -float("inf"):
            p = 10
        print("Walking Sum,cnt is {}--{}--Percent is {}---->P is {}".format( sumrow,len(data),percent,p ) )
        ds_pd = pd.DataFrame(data=data,columns=["Row","Column","Length"])
        ds_pd = ds_pd.sort_values("Length",ascending=False)
        tmp = ds_pd[ ds_pd.Length>limit ]
        for item in data:
            row,col,val = int(item[0]),int(item[1]),float(item[2])
            sumlen=val
            order,pq=[],[]
            if vis[row][col] == 0:
                vis[row][col] = 1
                Walking(row,col,ds_array,x,y,p,sumlen,limit )
        print("Route is ",len(route))
        tmp = pd.DataFrame(data = route,columns=["Route","SumLength"])
        tmp.to_csv("/home/xiatengxi/Documents/BigDataSystem_B/ReDeFineCity/"+name+"_ReDeFineCity.csv",index=False)
    endtime = datetime.datetime.now()
    print( (endtime - starttime).seconds,"(s) Finish")