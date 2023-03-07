# MIT License
#
# Copyright (c) 2020 Debopam Bhattacherjee
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import math
import ephem
import pandas as pd

try:
    from . import util
except (ImportError, SystemError):
    import util

# Visualizes paths between endpoints at specific time instances

EARTH_RADIUS = 6378135.0 # WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html

# CONSTELLATION GENERATION GENERAL CONSTANTS
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EPOCH = "2000-01-01 00:00:00"

# CONSTELLATION SPECIFIC PARAMETERS
# STARLINK 550
NAME = "starlink_550"

################################################################
# The below constants are taken from Starlink's FCC filing as below:
# [1]: https://fcc.report/IBFS/SAT-MOD-20190830-00087
################################################################

MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
ALTITUDE_M = 550000  # Altitude ~550 km
SATELLITE_CONE_RADIUS_M = 940700 # From https://fcc.report/IBFS/SAT-MOD-20181108-00083/1569860.pdf (minimum angle of elevation: 25 deg)
MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2)) # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
NUM_ORBS = 72
NUM_SATS_PER_ORB = 22
INCLINATION_DEGREE = 53


# KUIPER 630
"""
NAME = "kuiper_630"

################################################################
# The below constants are taken from Kuiper's FCC filing as below:
# [1]: https://www.itu.int/ITU-R/space/asreceived/Publication/DisplayPublication/8716
################################################################

MEAN_MOTION_REV_PER_DAY = 14.80  # Altitude ~630 km
ALTITUDE_M = 630000  # Altitude ~630 km
SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(30.0))  # Considering an elevation angle of 30 degrees; possible values [1]: 20(min)/30/35/45
MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))  # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
NUM_ORBS = 34
NUM_SATS_PER_ORB = 34
INCLINATION_DEGREE = 51.9
"""

# General files needed to generate visualizations; Do not change for different simulations
topFile = "../static_html/top.html"
bottomFile = "../static_html/bottom.html"
city_detail_file = "../../gen_data2/starlink_info/ground_stations_cities_sorted_by_estimated_2025_pop_top_1000.basic.txt"

new_york_2_los_path_file = "../../gen_data2/starlink_info/final_path_new_york_2_los.txt"

# Time in ms for which visualization will be generated
GEN_TIME=1  #ms

# Input file; Generated during simulation
# Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# City IDs are available in the city_detail_file.
# If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# then offset ID is 1584 + 24 = 1608.
path_file = "../../gen_data2/starlink_info/mesh_net_path.txt"
sats_nodes_file ='../../gen_data2/starlink_info/ground_station_satellites_in_range.txt'
# Output directory for creating visualization html files
OUT_DIR = "../viz_output/"
OUT_HTML_FILE = OUT_DIR + NAME + "_path"

sat_objs = []
city_details = {}
paths_over_time = []


def generate_path_at_time():
    """
    Generates end-to-end path at specified time
    :return: HTML formatted string for visualization
    """
    viz_string = ""
    global src_GS
    global dst_GS
    global paths_over_time
    global OUT_HTML_FILE

    with open(path_file, 'r') as file:
        line = file.readline()
        line = line.split(',')
        line.pop()
        mesh_node = [int(val) for val in line]

    shifted_epoch = (pd.to_datetime(EPOCH))
    print(shifted_epoch)

    in_range_sats = []
    two_sat = []
    with open(sats_nodes_file, 'r') as file:
        line = file.readline()
        line = line.split(',')[:-1]

        in_range_sats = [int(val) for val in line]
        # two_sat.append(in_range_sats[-2])
        # two_sat.append(in_range_sats[-1])
        # in_range_sats.pop()
        # in_range_sats.pop()

    path_new_york_2_los = []
    with open(new_york_2_los_path_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split(',')
            line = [int(val) for val in line]
            path_new_york_2_los.append(line)

    first_two_sataion = [path_new_york_2_los[0][0], path_new_york_2_los[0][-1]]
    second_two_sataion = [path_new_york_2_los[1][0], path_new_york_2_los[1][-1]]

    for i in range(len(sat_objs)):
        sat_objs[i]["sat_obj"].compute(shifted_epoch)
        if i in first_two_sataion:
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", " + str(sat_objs[i]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.BLACK.withAlpha(1),}});\n"
        elif i in second_two_sataion:
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", " + str(sat_objs[i]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.YELLOW.withAlpha(1),}});\n"
        elif i in in_range_sats:
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", " + str(
                sat_objs[i]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.PINK.withAlpha(1),}});\n"
        elif i in mesh_node:
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", " + str(sat_objs[i]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.RED.withAlpha(1),}});\n"
        else:
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", "+str(sat_objs[i]["alt_km"]*1000)+"), "\
                     + "ellipsoid : {radii : new Cesium.Cartesian3(20000.0, 20000.0, 20000.0), "\
                     + "material : Cesium.Color.BLUE.withAlpha(1),}});\n"


    orbit_links = util.find_orbit_links(sat_objs, NUM_ORBS, NUM_SATS_PER_ORB)
    for key in orbit_links:
        sat1 = orbit_links[key]["sat1"]
        sat2 = orbit_links[key]["sat2"]
        viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                      + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                      + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                      + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                      + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                      + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                      + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                      + "width: 0.5, arcType: Cesium.ArcType.NONE, " \
                      + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                      + "color: Cesium.Color.GREY.withAlpha(0.3), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    OUT_HTML_FILE += "-end-to-end"

    # 9,New-York-Newark,40.717042,-74.003663,0
    # 20,Los-Angeles-Long-Beach-Santa-Ana,34.031656,-118.241716,0
    two_citys = [9, 20]

    # 123, San - Francisco - Oakland, 37.759881, -122.437392, 0
    # 32, Bogotá, 4.60971, -74.08175, 0
    # two_citys = [123, 32]





    for city_id in two_citys:
        viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                              + str(city_details[city_id]["long_deg"]) + ", " \
                              + str(city_details[city_id]["lat_deg"]) + ", " \
                              + str(city_details[city_id]["alt_km"] * 1000) + "), " \
                              + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                              + "material : Cesium.Color.GREEN.withAlpha(1),}});\n"
    # sat 0
    src = path_new_york_2_los[0][0]
    dst = path_new_york_2_los[0][-1]
    viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                  + str(city_details[two_citys[0]]["long_deg"]) + "," \
                  + str(city_details[two_citys[0]]["lat_deg"]) + "," \
                  + str(city_details[two_citys[0]]["alt_km"] * 1000) + "," \
                  + str(math.degrees(sat_objs[src]["sat_obj"].sublong)) + "," \
                  + str(math.degrees(sat_objs[src]["sat_obj"].sublat)) + "," \
                  + str(sat_objs[src]["alt_km"] * 1000) + "]), " \
                  + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                  + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                  + "color: Cesium.Color.BLUE.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                  + str(city_details[two_citys[1]]["long_deg"]) + "," \
                  + str(city_details[two_citys[1]]["lat_deg"]) + "," \
                  + str(city_details[two_citys[1]]["alt_km"] * 1000) + "," \
                  + str(math.degrees(sat_objs[dst]["sat_obj"].sublong)) + "," \
                  + str(math.degrees(sat_objs[dst]["sat_obj"].sublat)) + "," \
                  + str(sat_objs[dst]["alt_km"] * 1000) + "]), " \
                  + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                  + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                  + "color: Cesium.Color.BLUE.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
    # sat1
    src1 = path_new_york_2_los[1][0]
    dst1 = path_new_york_2_los[1][-1]
    viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                  + str(city_details[two_citys[0]]["long_deg"]) + "," \
                  + str(city_details[two_citys[0]]["lat_deg"]) + "," \
                  + str(city_details[two_citys[0]]["alt_km"] * 1000) + "," \
                  + str(math.degrees(sat_objs[src1]["sat_obj"].sublong)) + "," \
                  + str(math.degrees(sat_objs[src1]["sat_obj"].sublat)) + "," \
                  + str(sat_objs[src1]["alt_km"] * 1000) + "]), " \
                  + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                  + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                  + "color: Cesium.Color.BLUE.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                  + str(city_details[two_citys[1]]["long_deg"]) + "," \
                  + str(city_details[two_citys[1]]["lat_deg"]) + "," \
                  + str(city_details[two_citys[1]]["alt_km"] * 1000) + "," \
                  + str(math.degrees(sat_objs[dst1]["sat_obj"].sublong)) + "," \
                  + str(math.degrees(sat_objs[dst1]["sat_obj"].sublat)) + "," \
                  + str(sat_objs[dst1]["alt_km"] * 1000) + "]), " \
                  + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                  + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                  + "color: Cesium.Color.BLUE.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    color_line_four_path = ['RED', 'ORANGE', 'BROWN', 'YELLOW']
    for idx, path in enumerate(path_new_york_2_los):
        for id in range(len(path) - 1):
            src = path[id]
            dst = path[id + 1]
            viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublat)) + "," + str(sat_objs[src]["alt_km"] * 1000) + "," \
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublat)) + "," + str(sat_objs[dst]["alt_km"] * 1000) + "]), " \
                          + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                          + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                          + "color: Cesium.Color.{}".format(color_line_four_path[idx]) + ".withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"


    OUT_HTML_FILE += "_" + str(GEN_TIME) + ".html"
    print('output_file', OUT_HTML_FILE)
    return viz_string


city_details = util.read_city_details(city_details, city_detail_file)
sat_objs = util.generate_sat_obj_list(
    NUM_ORBS,
    NUM_SATS_PER_ORB,
    EPOCH,
    PHASE_DIFF,
    INCLINATION_DEGREE,
    ECCENTRICITY,
    ARG_OF_PERIGEE_DEGREE,
    MEAN_MOTION_REV_PER_DAY,
    ALTITUDE_M
)
viz_string = generate_path_at_time()
util.write_viz_files(viz_string, topFile, bottomFile, OUT_HTML_FILE)
