import sys
sys.path.append("../satgenpy")
import satgen
import os
import networkx as nx
from astropy import units as u
import math
import shutil
from loguru import logger
import random


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

if os.path.exists(gen_data):
    shutil.rmtree(gen_data)
os.mkdir(gen_data)
# gen tle data
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

isl_file_all = os.path.join(gen_data, "isls")
if os.path.exists(isl_file_all):
    shutil.rmtree(isl_file_all)
os.mkdir(isl_file_all)

distance_file_all = os.path.join(gen_data, "distance")
if os.path.exists(distance_file_all):
    shutil.rmtree(distance_file_all)
os.mkdir(distance_file_all)

access_sat_file_all = os.path.join(gen_data, "access_sat")

if os.path.exists(access_sat_file_all):
    shutil.rmtree(access_sat_file_all)
    # os.mkdir(access_sat_file_all)
os.mkdir(access_sat_file_all)

sat_info = satgen.read_tles(os.path.join(gen_data, "tles.txt"))
satellites = sat_info['satellites']
epoch = sat_info['epoch']
isl_list = satgen.read_isls(os.path.join(gen_data, "isls.txt"), len(satellites))

ground_station_satellites_in_range = []

satgen.extend_ground_stations(
    "input_data/wlh.txt",
    gen_data + "/extended_ground_stations.txt")

print('start calculating ')
ground_stations = satgen.read_ground_stations_extended(gen_data + "/extended_ground_stations.txt")
for each_step in range(0, simulation_end_time_s, Time_step_s):
    each_time_current = epoch + each_step * u.second

    # write all sats distance to distance_.txt
    distance_file = os.path.join(distance_file_all, "distance_{}.txt".format(each_step))
    eachtime_isls_file = os.path.join(isl_file_all, "isls_{}.txt".format(each_step))
    with open(eachtime_isls_file, 'w') as f :
        f.write('src des bd delay plr \n')
    # with open(distance_file, 'w') as f:
    #     for (a, b) in isl_list:
    #         sat_distance_m = satgen.distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(each_time_current))
    #         f.write(str(sat_distance_m))
    #         f.write('\n')

    # calculate the shortest distance left and right sats
    isls_list = []
    for i in range(len(satellites)):
        # 当前卫星所属的轨道编号，和轨道内编号
        num_of_plane = i//NUM_SATS_PER_ORB
        num_of_in_plane = i % NUM_SATS_PER_ORB

        # 左边轨道相连接的最短距离的卫星 计算一圈
        left_plane = num_of_plane-1 if num_of_plane >= 1 else NUM_ORBS-1
        right_plane = (num_of_plane+1) % NUM_ORBS
        num_of_left_sat = left_plane * NUM_SATS_PER_ORB + num_of_in_plane
        num_of_right_sat = right_plane * NUM_SATS_PER_ORB + num_of_in_plane
        min_dist_left = satgen.distance_m_between_satellites(satellites[i], satellites[num_of_left_sat], str(epoch), str(each_time_current))
        min_dist_right = satgen.distance_m_between_satellites(satellites[i], satellites[num_of_right_sat], str(epoch), str(each_time_current))

        for sat_of_left in range(left_plane * NUM_SATS_PER_ORB, left_plane * NUM_SATS_PER_ORB + NUM_SATS_PER_ORB):
            sat_distance_m = satgen.distance_m_between_satellites(satellites[i], satellites[sat_of_left], str(epoch), str(each_time_current))
            if sat_distance_m < min_dist_left:
                min_dist_left = sat_distance_m
                num_of_left_sat = sat_of_left
        for sat_of_right in range(right_plane * NUM_SATS_PER_ORB, right_plane * NUM_SATS_PER_ORB + NUM_SATS_PER_ORB):
            sat_distance_m = satgen.distance_m_between_satellites(satellites[i], satellites[sat_of_right], str(epoch), str(each_time_current))
            if sat_distance_m < min_dist_right:
                min_dist_right = sat_distance_m
                num_of_right_sat = sat_of_right

        # 写isl连接关系到文件中 isls_xxx.txt 当前卫星的上下左右
        num_of_up_sat = num_of_plane * NUM_SATS_PER_ORB + num_of_in_plane-1 if num_of_in_plane >= 1 else NUM_SATS_PER_ORB-1
        num_of_down_sat = num_of_plane * NUM_SATS_PER_ORB + ((num_of_in_plane+1) % NUM_SATS_PER_ORB)

        with open(eachtime_isls_file, 'a') as f:
            if (i, num_of_up_sat) not in isls_list and (num_of_up_sat, i) not in isls_list:
                isls_list.append((i, num_of_up_sat))
                bd = round(random.uniform(100, 110), 6)
                delay = round(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_up_sat], str(epoch),
                                                                   str(each_time_current)) / 3e8, 6)
                plr = round(random.uniform(0.005, 0.015), 6)
                f.write('{} {} {} {} {}'.format(i, num_of_up_sat, bd, delay, plr))
                f.write('\n')

            if (i, num_of_down_sat) not in isls_list and (num_of_down_sat, i) not in isls_list:
                isls_list.append((i, num_of_down_sat))
                bd = round(random.uniform(100, 110), 6)
                delay = round(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_down_sat], str(epoch),
                                                                   str(each_time_current)) / 3e8, 6)
                plr = round(random.uniform(0.005, 0.015), 6)
                f.write('{} {} {} {} {}'.format(i, num_of_down_sat, bd, delay, plr))
                f.write('\n')

            if (i, num_of_left_sat) not in isls_list and (num_of_left_sat, i) not in isls_list:
                isls_list.append((i, num_of_left_sat))
                bd = round(random.uniform(100, 110), 6)
                delay = round(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_left_sat], str(epoch),
                                                                   str(each_time_current)) / 3e8, 6)
                plr = round(random.uniform(0.005, 0.015), 6)
                f.write('{} {} {} {} {}'.format(i, num_of_left_sat, bd, delay, plr))
                f.write('\n')

            if (i, num_of_right_sat) not in isls_list and (num_of_right_sat, i) not in isls_list:
                isls_list.append((i, num_of_right_sat))
                bd = round(random.uniform(100, 110), 6)
                delay = round(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_right_sat], str(epoch),
                                                                   str(each_time_current)) / 3e8, 6)
                plr = round(random.uniform(0.005, 0.015), 6)
                f.write('{} {} {} {} {}'.format(i, num_of_right_sat, bd, delay, plr))
                f.write('\n')
        # 写距离到文件 distance_{}.txt
        with open(distance_file, 'a') as f:
            f.write(str(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_up_sat], str(epoch), str(each_time_current))))
            f.write('\n')
            f.write(str(satgen.distance_m_between_satellites(satellites[i], satellites[num_of_down_sat], str(epoch), str(each_time_current))))
            f.write('\n')
            f.write(str(min_dist_left))
            f.write('\n')
            f.write(str(min_dist_right))
            f.write('\n')
    isls_list.clear()
    # 计算所有卫星的距离并保存到文件中
    # with open(distance_file, 'w') as f:
    #     for i in range(72*22):
    #         for j in range(72*22):
    #             if i == j:
    #                 sat_distance_m = 0
    #             else:
    #                 sat_distance_m = satgen.distance_m_between_satellites(satellites[i], satellites[j], str(epoch),
    #                                                                       str(each_time_current))
    #             f.write(str(sat_distance_m))
    #             f.write(',')
    #         f.write('\n')

    # write access sat to ground
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
                # if distance_m <= MAX_GSL_LENGTH_M:
                satellites_in_range.append((distance_m, sid))

            satellites_in_range = sorted(satellites_in_range)
            # logger.info(satellites_in_range)
            f.write(str(satellites_in_range))
            f.write('\n')


logger.info("finish")