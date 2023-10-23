import sys
sys.path.append("../satgenpy")
import satgen
import os
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


sat_info = satgen.read_tles(os.path.join(gen_data, "tles.txt"))
satellites = sat_info['satellites']
epoch = sat_info['epoch']


def get_sat_lon_and_lat(sats, epoch):
    for idx, sat in enumerate(sats):
        sat.compute(str(epoch))
        satellite_longitude = sat.sublong
        satellite_latitude = sat.sublat
        with open(os.path.join(gen_data, 'drgee.txt'), 'a') as f:
            f.write(str(satellite_latitude))
            f.write('\n')
            f.write(str(satellite_longitude))
            f.write('\n')

get_sat_lon_and_lat(satellites, epoch)