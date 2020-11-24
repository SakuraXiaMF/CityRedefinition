import numpy as np

from osgeo import ogr
import osr
import gdal
from shapely.geometry import Point, shape
import geopandas as gpd

import os
import math
import multiprocessing


def OSMtoShape(osm_bz_path, shape_path=None, query=None):
    os.system("bzcat {0} | ogr2ogr -f SQLite _temp.sqlite /vsistdin/".format(osm_bz_path))


def mergeShapefile(shp_list, merged_shp):
    suffix = [".dbf", ".prj", ".shp", ".shx"]

    for f in shp_list:
        if os.path.exists(merged_shp):
            os.system("ogr2ogr -f 'ESRI Shapefile' -update -append {0} {1}".format(merged_shp, f))
        else:
            f_dir = os.path.dirname(f)
            f_base = os.path.splitext(os.path.basename(f))[0]
            merged_shp_dir = os.path.dirname(merged_shp)
            merged_shp_base = os.path.splitext(os.path.basename(merged_shp))[0]
            for s in suffix:
                os.system("cp {0} {1}".format(os.path.join(f_dir, f_base + s), os.path.join(merged_shp_dir, merged_shp_base + s)))


def createFishNet(input_shp, output_grid, height, width, driver_name="ESRI Shapefile", extent=None):
    # ------open and read Shapefile
    driver = ogr.GetDriverByName(driver_name)
    shp_ds = driver.Open(input_shp, 1)
    shp_layer = shp_ds.GetLayer()

    # ------get features' extent (or use specified extent) and projection information
    if extent is None:
        x_min, x_max, y_min, y_max = shp_layer.GetExtent()
    else:
        x_min, y_min, x_max, y_max = extent

    input_proj = shp_layer.GetSpatialRef()

    # ------define x,y coordinates of output FishNet
    num_row = math.ceil((y_max - y_min) / height)
    num_col = math.ceil((x_max - x_min) / width)

    fishnet_X_left = np.linspace(x_min, x_min+(num_col-1)*width, num_col)
    fishnet_X_right = fishnet_X_left + width
    fishnet_Y_top = np.linspace(y_max-(num_row-1)*width, y_max, num_row)
    fishnet_Y_top = np.ascontiguousarray(fishnet_Y_top[::-1])
    fishnet_Y_bottom = fishnet_Y_top - height

    # ------create output file
    out_driver = ogr.GetDriverByName(driver_name)
    if os.path.exists(output_grid):
        os.remove(output_grid)
    out_ds = out_driver.CreateDataSource(output_grid)
    out_layer = out_ds.CreateLayer(output_grid, geom_type=ogr.wkbPolygon)
    feature_def = out_layer.GetLayerDefn()

    # ------create features of grid cell
    for i in range(0, num_row):
        y_top = fishnet_Y_top[i]
        y_bottom = fishnet_Y_bottom[i]
        for j in range(0, num_col):
            x_left = fishnet_X_left[j]
            x_right = fishnet_X_right[j]
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x_left, y_top)
            ring.AddPoint(x_right, y_top)
            ring.AddPoint(x_right, y_bottom)
            ring.AddPoint(x_left, y_bottom)
            ring.AddPoint(x_left, y_top)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            #poly.Transform(coord_trans)
            # ---------add new geometry to layer
            out_feature = ogr.Feature(feature_def)
            out_feature.SetGeometry(poly)
            out_layer.CreateFeature(out_feature)
            out_feature.Destroy()

    # ------close data sources
    out_ds.Destroy()

    # ------write ESRI.prj file
    input_proj.MorphToESRI()
    file = open(os.path.splitext(output_grid)[0]+".prj", 'w')
    file.write(input_proj.ExportToWkt())
    file.close()


# ************************* [1] Road Network Intersection Points *************************
def createRoadIntersectionPoint(road_shp, output_shp, suffix=None):
    if suffix is not None:
        output_dir = os.path.dirname(output_shp)
        output_base = os.path.splitext(os.path.basename(road_shp))[0]
        output_shp = os.path.join(output_dir, output_base + suffix + ".shp")

    road_ds = ogr.Open(road_shp)
    road_layer = road_ds.GetLayer()
    x_min, x_max, y_min, y_max = road_layer.GetExtent()

    dy = y_max - y_min
    dx = x_max - x_min

    if dy * dx != 0:
        res_list = []
        roads = gpd.read_file(road_shp)
        for r_id, r_row in roads.iterrows():
            line_other = [shape(row["geometry"]) for index, row in roads.iterrows() if index > r_id]
            gdf_other = gpd.GeoDataFrame(geometry=line_other)
            gdf_r = gpd.GeoDataFrame(geometry=[r_row["geometry"]])
            intersect_pt = gdf_r.unary_union.intersection(gdf_other.unary_union)

            # ---ref to: https://gis.stackexchange.com/questions/137909/intersecting-lines-to-get-crossings-using-python-with-qgis
            if intersect_pt.type == "Point":
                res_list.append(intersect_pt)
            elif intersect_pt.type == "MultiPoint":
                res_list.extend([pt for pt in intersect_pt])
            elif intersect_pt.type == "LineString":
                if not intersect_pt.is_empty:
                    res_list.append(Point(intersect_pt.coords[0]))
                    res_list.append(Point(intersect_pt.coords[-1]))
            elif intersect_pt.type == "LineString":
                for geom in intersect_pt:
                    if geom.type == "Point":
                        res_list.append(intersect_pt)
                    elif geom.type == "MultiPoint":
                        res_list.extend([pt for pt in intersect_pt])
                    elif geom.type == "LineString":
                        if not intersect_pt.is_empty:
                            res_list.append(Point(intersect_pt.coords[0]))
                            res_list.append(Point(intersect_pt.coords[-1]))

        gdf_res = gpd.GeoDataFrame(geometry=res_list)
        gdf_res.to_file(output_shp)
        res = output_shp
    else:
        res = "EMPTY"

    return res


def createRoadIntersectionPoint_multiproc(road_shp, output_shp, reserved=False, num_cpu=1):
    # ------get the basic information of the road network layer
    road_dir = os.path.dirname(road_shp)
    road_base = os.path.splitext(os.path.basename(road_shp))[0]
    road_ds = ogr.Open(road_shp)
    road_layer = road_ds.GetLayer()
    x_min, x_max, y_min, y_max = road_layer.GetExtent()

    # ------create sub-regions for parallel processing
    num_cpu_available = multiprocessing.cpu_count()
    if num_cpu_available < num_cpu:
        print("Use %d CPUs which is available now" % num_cpu_available)
        num_cpu = num_cpu_available

    n_sub = math.floor(math.sqrt(num_cpu))
    dx = (x_max - x_min) / n_sub
    dy = (y_max - y_min) / n_sub
    x_min_list = np.array([x_min + i * dx for i in range(0, n_sub)])
    x_max_list = x_min_list + dx
    y_min_list = np.array([y_min + i * dy for i in range(0, n_sub)])
    y_max_list = y_min_list + dy

    arg_list = []
    for i in range(0, n_sub):
        for j in range(0, n_sub):
            subRegion = os.path.join(road_dir, road_base + "_temp_{0}.shp".format(str(i) + str(j)))
            subOutput = os.path.join(road_dir, road_base + "_IntersectionPt_temp_{0}.shp".format(str(i) + str(j)))
            extent = [x_min_list[j], y_min_list[i], x_max_list[j], y_max_list[i]]
            os.system("ogr2ogr -f 'ESRI Shapefile' {0} {1} -clipsrc {2}".format(subRegion, road_shp, " ".join(str(x) for x in extent)))
            arg_list.append((subRegion, subOutput))

    # ------call the createRoadIntersectionPoint function for each sub-regions
    pool = multiprocessing.Pool(processes=num_cpu)
    res_list = pool.starmap(createRoadIntersectionPoint, arg_list)
    pool.close()
    pool.join()

    # ------merge GeoTiff from sub-regions
    res_shp_list = [i for i in res_list if i != "EMPTY"]
    mergeShapefile(shp_list=res_shp_list, merged_shp=output_shp)

    if not reserved:
        os.system("rm {0}".format(os.path.join(road_dir, "*_temp_*")))

    print("The Intersection Point Map has been created at: {0}".format(output_shp))


# ************************* [2] Road Network Length *************************
def getRoadLength(road_shp, output_grid, resolution=0.5/60.0, scale=1.0, reserved=False, suffix=None, extent=None):
    road_dir = os.path.dirname(road_shp)
    fishnet_name = os.path.join(road_dir, "_fishnet_{0}.shp".format(suffix))

    road_ds = ogr.Open(road_shp)
    road_layer = road_ds.GetLayer()
    x_min, x_max, y_min, y_max = road_layer.GetExtent()
    input_proj = road_layer.GetSpatialRef()

    num_row = math.ceil((y_max - y_min) / resolution)
    num_col = math.ceil((x_max - x_min) / resolution)
    print( "GetRoadLength----DataLength is " , len(output_grid) )
    if num_col * num_row != 0:
        # ------create FishNet layer for the entire Shapefile
        createFishNet(input_shp=road_shp, output_grid=fishnet_name, height=resolution, width=resolution, extent=extent)

        # ------get the intersection part of road network layer and FishNet layer
        # ------the output layer contains segmented roads using a field named "FN_FID" marking which cell each part belongs to
        fishnet_ds = ogr.Open(fishnet_name)
        fishnet_layer = fishnet_ds.GetLayer()

        intersect_driver = ogr.GetDriverByName("ESRI Shapefile")
        intersect_path = os.path.join(road_dir, "_intersect_{0}.shp".format(suffix))
        if os.path.exists(intersect_path):
            os.remove(intersect_path)
        print("Create Intersect_DS")
        intersect_ds = intersect_driver.CreateDataSource(intersect_path)
        intersect_layer = intersect_ds.CreateLayer(output_grid, geom_type=ogr.wkbLineString)
        print("road_Layer Intersection()")
        road_layer.Intersection(fishnet_layer, intersect_layer, ["METHOD_PREFIX=FN_"])

        # ------calculate the sum of road length for each cell in FishNet layer
        val_list = [[feature.GetField("FN_FID"), feature.GetGeometryRef().Length()] for feature in intersect_layer]

        if extent is not None:
            x_min, y_min, x_max, y_max = extent
            num_row = math.ceil((y_max - y_min) / resolution)
            num_col = math.ceil((x_max - x_min) / resolution)

        length_arr = np.zeros(num_row * num_col)
        for v in val_list:
            length_arr[v[0]] += v[1] * scale

        length_arr = length_arr.reshape((num_row, num_col))

        driver = gdal.GetDriverByName("GTiff")
        print("Driver.Create")
        length_ds = driver.Create(output_grid, num_col, num_row, 1, gdal.GDT_Float64)
        print("GetBand")
        length_ds.GetRasterBand(1).WriteArray(length_arr)
        print("SetSomething")
        length_ds.SetGeoTransform([x_min, resolution, 0, y_max, 0, -resolution])
        length_ds.SetProjection(input_proj.ExportToWkt())
        print("FlushCache")
        length_ds.FlushCache()
        length_ds = None
        # ---------if reserve_flag is set to be True, we store intermediate results for further algorithm validation
        if not reserved:
            os.system("rm {0}".format(os.path.join(road_dir, "_*")))
        res = output_grid
    else:
        res = "EMPTY"

    return res


def getRoadLength_multiproc(road_shp, resolution=0.5/60.0, scale=1.0, reserved=False, num_cpu=1):
    # ------get the basic information of the road network layer
    road_dir = os.path.dirname(road_shp)
    road_base = os.path.splitext(os.path.basename(road_shp))[0]
    road_ds = ogr.Open(road_shp)
    road_layer = road_ds.GetLayer()
    x_min, x_max, y_min, y_max = road_layer.GetExtent()

    # ------create sub-regions for parallel processing
    num_cpu_available = multiprocessing.cpu_count()
    print("Use %d CPUs which is available now" % num_cpu_available)
    if num_cpu_available < num_cpu:
        print("Use %d CPUs which is available now" % num_cpu_available)
        num_cpu = num_cpu_available

    n_sub = math.floor(math.sqrt(num_cpu))
    dx = (x_max - x_min) / n_sub
    dy = (y_max - y_min) / n_sub
    x_min_list = np.array([x_min + i * dx for i in range(0, n_sub)])
    x_max_list = x_min_list + dx
    y_min_list = np.array([y_min + i * dy for i in range(0, n_sub)])
    y_max_list = y_min_list + dy

    arg_list = []
    print("{}Line {}Column Need to Finish".format(n_sub,n_sub))
    for i in range(0, n_sub): 
        print("************************** {} Line Dealing **************************".format(i))
        for j in range(0, n_sub):
            subRegion = os.path.join(road_dir, road_base + "_temp_{0}.shp".format(str(i) + str(j)))
            subOutput = os.path.join(road_dir, road_base + "_roadLength_temp_{0}.tif".format(str(i) + str(j)))
            extent = [x_min_list[j], y_min_list[i], x_max_list[j], y_max_list[i]]
            os.system("ogr2ogr -f 'ESRI Shapefile' {0} {1} -clipsrc {2}".format(subRegion, road_shp, " ".join(str(x) for x in extent)))
            arg_list.append((subRegion, subOutput, resolution, scale, True, str(i) + str(j), extent))
    
    print("************************** ESRI Shapefile Finish **************************")
    # ------call the conversion function for each sub-regions
    pool = multiprocessing.Pool(processes=num_cpu)
    res_list = pool.starmap(getRoadLength, arg_list)
    pool.close()
    pool.join()

    # ------merge GeoTiff from sub-regions
    print("road_dir",road_dir+"/TIF/")
    # print("road_base",road_base)
    tiff_path = os.path.join(road_dir+"/TIF/", road_base + "_roadLength.tif")

    res_tiff_list = [i for i in res_list if i != "EMPTY"]
    os.system("gdalwarp {0} {1}".format(" ".join(res_tiff_list), tiff_path))

    if not reserved:
        os.system("rm {0}".format(os.path.join(road_dir, "*_temp_*")))
        os.system("rm {0}".format(os.path.join(road_dir, "_*")))

    print("The Road Length Map has been created at: {0}".format(tiff_path))


from osgeo import gdal

def getsuffix(path,suffix):
    res = []
    datanames = os.listdir(path)
    for dataname in datanames:
        if os.path.splitext(dataname)[1] == suffix:
            res.append( path+dataname)
    return res


def NodealFiles():
    shpfiles =getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/",".shp")
    tiffiles =getsuffix("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/TIF/",".tif")
    SHP = "/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/IDN.shp"
    vis = {}
    for shp in tiffiles:
        name = os.path.basename(shp)[:-4]
        idx = name.find("_roadLength")
        vis[name[:idx]]=1

    shpidx = SHP.find(".shp")
    idx = shpfiles.index(SHP)
    # print( len(shpfiles),idx,os.path.getsize(shpfiles[idx+1])/1024/1024)

    NoDealFile=[]
    for shp in shpfiles:
        basename = os.path.basename(shp)
        name = os.path.basename(shp)[:-4]
        idx = name.find("_")
        if name == 'GEO' or name == "world" or name == "GHA":
            continue
        if idx==-1:
            if name not in vis:
                size = os.path.getsize(shp)/1024/1024/1024
                if size<1:
                    NoDealFile.append(shp)
                    # print(basename,end='\t')        
    # print('\n',len(NoDealFile),NoDealFile[0])
    return NoDealFile

if __name__ == "__main__":
    path = "/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/"
    # shapefile = getsuffix(path,".shp")
    # idx = shapefile.index("/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/DEU.shp")
    #Failes = ["/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/GEO.shp","/home/xiatengxi/Documents/BigDataSystem_B/WorldSHP/world.shp"]
    shpfiles = NodealFiles()
    f = open("./now.txt", "a")
    for shape in shpfiles:#shapefile[225:]:
        file = shape
        f.write(file+"\n")
        # f.close()
        print(file,"Start")
        road = os.path.join(file)
        getRoadLength_multiproc(road, reserved=False, num_cpu=15)
        print("Finish")
    f.close()