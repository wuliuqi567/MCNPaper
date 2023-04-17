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

# from astropy import units as u
# from poliastro.bodies import Earth
# from poliastro.twobody import Orbit
# from astropy.time import Time
# from extractor import CZMLExtractor
import math
try:
    from . import util
except (ImportError, SystemError):
    import util

# Generate static visualizations for entire constellation (multiple shells).

EARTH_RADIUS = 6378135.0 # WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html

# CONSTELLATION GENERATION GENERAL CONSTANTS
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EPOCH = "2000-01-01 00:00:00"

# Shell wise color codes
# COLOR = [[255, 0, 0, 200], [32, 128, 46, 200], [0, 0, 255, 200], [245, 66, 242, 200], [245, 126, 66, 200]]
COLOR = ['CRIMSON', 'FORESTGREEN', 'DODGERBLUE', 'PERU', 'BLUEVIOLET', 'DARKMAGENTA']
# CONSTELLATION SPECIFIC PARAMETERS


# STARLINK
NAME = "routing_mesh_net"

SHELL_CNTR = 5

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 15.19  # Altitude ~550000 km
ALTITUDE_M[0] = 550000  # Altitude ~550000 km
NUM_ORBS[0] = 72
NUM_SATS_PER_ORB[0] = 22
INCLINATION_DEGREE[0] = 53
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

MEAN_MOTION_REV_PER_DAY[1] = 13.4  # Altitude ~1110 km
ALTITUDE_M[1] = 1110000  # Altitude ~1110 km
NUM_ORBS[1] = 32
NUM_SATS_PER_ORB[1] = 50
INCLINATION_DEGREE[1] = 53.8
BASE_ID[1] = 1584
ORB_WISE_IDS[1] = []

MEAN_MOTION_REV_PER_DAY[2] = 13.35  # Altitude ~1130 km
ALTITUDE_M[2] = 1130000  # Altitude ~1130 km
NUM_ORBS[2] = 8
NUM_SATS_PER_ORB[2] = 50
INCLINATION_DEGREE[2] = 74
BASE_ID[2] = 3184
ORB_WISE_IDS[2] = []

MEAN_MOTION_REV_PER_DAY[3] = 12.97  # Altitude ~1275 km
ALTITUDE_M[3] = 1275000  # Altitude ~1275 km
NUM_ORBS[3] = 5
NUM_SATS_PER_ORB[3] = 75
INCLINATION_DEGREE[3] = 81
BASE_ID[3] = 3584
ORB_WISE_IDS[3] = []

MEAN_MOTION_REV_PER_DAY[4] = 12.84  # Altitude ~1325 km
ALTITUDE_M[4] = 1325000  # Altitude ~1325 km
NUM_ORBS[4] = 6
NUM_SATS_PER_ORB[4] = 75
INCLINATION_DEGREE[4] = 70
BASE_ID[4] = 3959
ORB_WISE_IDS[4] = []


"""
# TELESAT
NAME = "Telesat"
SHELL_CNTR = 2

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 13.66  # Altitude ~1015 km
ALTITUDE_M[0] = 1015000  # Altitude ~1015 km
NUM_ORBS[0] = 27
NUM_SATS_PER_ORB[0] = 13
INCLINATION_DEGREE[0] = 98.98
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

MEAN_MOTION_REV_PER_DAY[1] = 12.84  # Altitude ~1325 km
ALTITUDE_M[1] = 1325000  # Altitude ~1325 km
NUM_ORBS[1] = 40
NUM_SATS_PER_ORB[1] = 33
INCLINATION_DEGREE[1] = 50.88
BASE_ID[1] = 351
ORB_WISE_IDS[1] = []
"""

"""
# KUIPER
NAME = "kuiper"
################################################################
# The below constants are taken from Kuiper's FCC filing as below:
# [1]: https://www.itu.int/ITU-R/space/asreceived/Publication/DisplayPublication/8716
################################################################

SHELL_CNTR = 3

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 14.80  # Altitude ~630 km
ALTITUDE_M[0] = 630000  # Altitude ~630 km
NUM_ORBS[0] = 34
NUM_SATS_PER_ORB[0] = 34
INCLINATION_DEGREE[0] = 51.9
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

MEAN_MOTION_REV_PER_DAY[1] = 14.86  # Altitude ~610 km
ALTITUDE_M[1] = 610000  # Altitude ~610 km
NUM_ORBS[1] = 36
NUM_SATS_PER_ORB[1] = 36
INCLINATION_DEGREE[1] = 42
BASE_ID[1] = 1156
ORB_WISE_IDS[1] = []

MEAN_MOTION_REV_PER_DAY[2] = 14.93  # Altitude ~590 km
ALTITUDE_M[2] = 590000  # Altitude ~590 km
NUM_ORBS[2] = 28
NUM_SATS_PER_ORB[2] = 28
INCLINATION_DEGREE[2] = 33
BASE_ID[2] = 2452
ORB_WISE_IDS[2] = []
"""


# General files needed to generate visualizations; Do not change for different simulations
topFile = "../static_html/top.html"
bottomFile = "../static_html/bottom.html"

# Output directory for creating visualization html files
OUT_DIR = "../viz_output/"
# JSON_NAME  = NAME+"_5shell.json"
# OUT_JSON_FILE = OUT_DIR + JSON_NAME
OUT_HTML_FILE = OUT_DIR + NAME + ".html"

# START = Time(EPOCH, scale="tdb")
# END = START + (10*60) * u.second
# sample_points = 10
# extractor = CZMLExtractor(START, END, sample_points)


def generate_satellite_trajectories():
    """
    Generates and adds satellite orbits to visualization.
    :return: viz_string
    """
    viz_string = ""
    for i in range(0, 1):
        sat_objs = util.generate_sat_obj_list(
            NUM_ORBS[i],
            NUM_SATS_PER_ORB[i],
            EPOCH,
            PHASE_DIFF,
            INCLINATION_DEGREE[i],
            ECCENTRICITY,
            ARG_OF_PERIGEE_DEGREE,
            MEAN_MOTION_REV_PER_DAY[i],
            ALTITUDE_M[i]
        )

        sats_in_mesh_net = get_mesh_net(1521, 903)["sats_in_mesh_list"]
        # for j in range(len(sat_objs)):
        for j in sats_in_mesh_net:
            sat_objs[j]["sat_obj"].compute(EPOCH)
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[j]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[j]["sat_obj"].sublat)) + ", " + str(
                sat_objs[j]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(30000.0, 30000.0, 30000.0), " \
                          + "material : Cesium.Color.BLACK.withAlpha(1),}});\n"
        orbit_links = util.find_grid_links(sat_objs, NUM_ORBS[i], NUM_SATS_PER_ORB[i])

        for key in orbit_links:
            sat1 = orbit_links[key]["sat1"]
            sat2 = orbit_links[key]["sat2"]
            if sat1 in sats_in_mesh_net and sat2 in sats_in_mesh_net:
                viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                              + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                              + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                              + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                              + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                              + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                              + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                              + "width: 0.5, arcType: Cesium.ArcType.NONE, " \
                              + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                              + "color: Cesium.Color."+COLOR[i]+".withAlpha(0.4), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
    return viz_string


def ascending_or_descending_of_satellite(satellite, date_str) -> int:
    satellite.compute(date_str)
    meanAnomaly = satellite.M  # Mean anomaly from perigee at epoch
    if 90 < meanAnomaly <= 270:
        # descending satellite
        return 0
    else:
        # ascending satellite
        return 1

def get_mesh_net(sat_src, sat_des):
    n_orbits = 72
    n_sats_per_orbit = 22

    sat_two_final_selected = []
    sat_two_final_selected.append([sat_src // n_sats_per_orbit, sat_src % n_sats_per_orbit])
    sat_two_final_selected.append([sat_des // n_sats_per_orbit, sat_des % n_sats_per_orbit])


    sat_connect2gs0_id = sat_two_final_selected[0]
    sat_connect2gs1_id = sat_two_final_selected[1]
    # three situation to specifly direction
    flag = 0
    if sat_connect2gs0_id[0] == sat_connect2gs1_id[0]:
        flag = 1
    elif sat_connect2gs0_id[1] == sat_connect2gs1_id[1]:
        flag = 2

    # bulid mesh network
    orbit_begin = sat_connect2gs0_id[0]
    orbit_end = sat_connect2gs1_id[0]

    if flag != 2:
        if orbit_begin < orbit_end:
            orbit_begin = (orbit_begin - 1) % n_orbits
            orbit_end = (orbit_end + 1) % n_orbits
        elif orbit_begin > orbit_end:
            orbit_end = (orbit_end - 1) % n_orbits
            orbit_begin = (orbit_begin + 1) % n_orbits
            orbit_begin, orbit_end = orbit_end, orbit_begin
        else:
            orbit_begin = (orbit_begin - 1) % n_orbits
            orbit_end = (orbit_end + 1) % n_orbits
    else:
        orbit_begin = min(sat_connect2gs0_id[0], sat_connect2gs1_id[0])
        orbit_end = max(sat_connect2gs0_id[0], sat_connect2gs1_id[0])

    n_sat_of_plane_begin = sat_connect2gs0_id[1]
    n_sat_of_plane_end = sat_connect2gs1_id[1]

    if flag != 1:
        if n_sat_of_plane_begin < n_sat_of_plane_end:
            n_sat_of_plane_begin = (n_sat_of_plane_begin - 1) % n_sats_per_orbit
            n_sat_of_plane_end = (n_sat_of_plane_end + 1) % n_sats_per_orbit
        elif n_sat_of_plane_begin > n_sat_of_plane_end:
            n_sat_of_plane_begin = (n_sat_of_plane_begin + 1) % n_sats_per_orbit
            n_sat_of_plane_end = (n_sat_of_plane_end - 1) % n_sats_per_orbit
            n_sat_of_plane_begin, n_sat_of_plane_end = n_sat_of_plane_end, n_sat_of_plane_begin
        else:
            n_sat_of_plane_begin = (n_sat_of_plane_begin - 1) % n_sats_per_orbit
            n_sat_of_plane_end = (n_sat_of_plane_end + 1) % n_sats_per_orbit
    else:
        n_sat_of_plane_begin = min(sat_connect2gs0_id[1], sat_connect2gs1_id[1])
        n_sat_of_plane_end = max(sat_connect2gs0_id[1], sat_connect2gs1_id[1])

    orbits_list = [orbit_begin]
    n_sat_plane_list = [n_sat_of_plane_begin]

    temp_orbit = ((orbit_begin + 1) % n_orbits)
    while (temp_orbit != orbit_end):
        orbits_list.append(temp_orbit)
        temp_orbit = ((temp_orbit + 1) % n_orbits)
    orbits_list.append(orbit_end)

    temp_sat_num = (n_sat_of_plane_begin + 1) % n_sats_per_orbit
    while (temp_sat_num != n_sat_of_plane_end):
        n_sat_plane_list.append(temp_sat_num)
        temp_sat_num = ((temp_sat_num + 1) % n_sats_per_orbit)
    n_sat_plane_list.append(n_sat_of_plane_end)

    # mesh_net = nx.Graph()
    sats_in_mesh_list = []
    for x in orbits_list:
        for y in n_sat_plane_list:
            sat_id = (x * n_sats_per_orbit + y)
            # mesh_net.add_node(sat_id, pos=(x, y))
            sats_in_mesh_list.append(sat_id)
    net_info = {
        "orbits_list": orbits_list,
        "n_sat_plane_list": n_sat_plane_list,
        "sats_in_mesh_list": sats_in_mesh_list
    }

    return net_info


def write_viz_files():
    """
    Writes JSON and TML files to the output folder
    :return: None
    """
    writer_html = open(OUT_HTML_FILE, 'w')
    with open(topFile, 'r') as fi:
        writer_html.write(fi.read())
    writer_html.write(viz_string)
    with open(bottomFile, 'r') as fb:
        writer_html.write(fb.read())
    writer_html.close()


viz_string = generate_satellite_trajectories()

write_viz_files()