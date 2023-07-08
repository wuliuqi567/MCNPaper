import sys
sys.path.append("./satgenpy")
import satgen
import os
import networkx as nx
from astropy import units as u
import math
import shutil
from loguru import logger
# read configure file
configure_dict = {}

with open('config.txt', 'r') as config_file:
    lines = config_file.readlines()
    for line in lines:
        line = line.strip().split('=')
        configure_dict.setdefault(line[0], line[1])

gen_data = configure_dict['satellite_network_dir_and_name']
ALTITUDE_M = int(configure_dict['ALTITUDE_M'])
NUM_ORBS = int(configure_dict['NUM_ORBS'])
NUM_SATS_PER_ORB = int(configure_dict['NUM_SATS_PER_ORB'])
INCLINATION_DEGREE = int(configure_dict['INCLINATION_DEGREE'])
simulation_end_time_s = int(configure_dict['simulation_end_time_s'])
Time_step_s = int(configure_dict['Time_step_s'])
Elevation_angle = int(configure_dict['Elevation_angle'])
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EARTH_RADIUS = 6378135.0
T = pow((ALTITUDE_M + EARTH_RADIUS) / 21613.546, 3 / 2)
MEAN_MOTION_REV_PER_DAY = round(24 * 60 * 60 / T, 2)
# Considering an elevation angle of 30 degrees; possible values [1]: 20(min)/30/35/45
SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(Elevation_angle))

MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))

if not os.path.exists(gen_data):
    os.makedirs(gen_data)

satgen.generate_tles_from_scratch_manual(
            gen_data + "/tles.txt",
            gen_data,
            NUM_ORBS,
            NUM_SATS_PER_ORB,
            PHASE_DIFF,
            INCLINATION_DEGREE,
            ECCENTRICITY,
            ARG_OF_PERIGEE_DEGREE,
            MEAN_MOTION_REV_PER_DAY
        )

        # ISLs
print("Generating ISLs...")

satgen.generate_plus_grid_isls(
    gen_data + "/isls.txt",
    NUM_ORBS,
    NUM_SATS_PER_ORB,
    isl_shift=0,
    idx_offset=0
)

tle_file = os.path.join(gen_data, "tles.txt")
isl_file = os.path.join(gen_data, "isls.txt")
distance_file_all = os.path.join(gen_data, "distance")
if os.path.exists(distance_file_all):
    shutil.rmtree(distance_file_all)
    os.mkdir(distance_file_all)

access_sat_file_all = os.path.join(gen_data, "access_sat")

print(access_sat_file_all)
if not os.path.exists(access_sat_file_all):
    # shutil.rmtree(access_sat_file_all)
    os.mkdir(access_sat_file_all)

sat_info = satgen.read_tles(tle_file)
satellites = sat_info['satellites']
epoch = sat_info['epoch']

isl_list = satgen.read_isls(isl_file, len(satellites))

ground_station_satellites_in_range = []

satgen.extend_ground_stations(
    "input_data/wlh.txt",
    gen_data + "/extended_ground_stations.txt")

ground_stations = satgen.read_ground_stations_extended("dsj_starlink/extended_ground_stations.txt")
for each_step in range(0, simulation_end_time_s, Time_step_s):
    each_time_current = epoch + each_step * u.second
    distance_file = os.path.join(distance_file_all, "distance_{}.txt".format(each_step))
    with open(distance_file, 'w') as f:
        for (a, b) in isl_list:
            sat_distance_m = satgen.distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(each_time_current))
            f.write(str(sat_distance_m))
            f.write('\n')

    access_sat_file = os.path.join(access_sat_file_all, "access_sat_{}.txt".format(each_step))
    with open(access_sat_file, 'w') as f:
        for ground_station in ground_stations:
            # Find satellites in range
            satellites_in_range = []
            for sid in range(len(satellites)):
                distance_m = satgen.distance_m_ground_station_to_satellite(
                    ground_station,
                    satellites[sid],
                    str(epoch),
                    str(each_time_current)
                )
                if distance_m <= MAX_GSL_LENGTH_M:
                    satellites_in_range.append((distance_m, sid))

            satellites_in_range = sorted(satellites_in_range)
            # logger.info(satellites_in_range)
            f.write(str(satellites_in_range))
            f.write('\n')


logger.info("finish")