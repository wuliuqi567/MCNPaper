import sys
sys.path.append("./satgenpy")
import satgen
import os
import networkx as nx
# satellite_network_dir_and_name=globalstar
# ALTITUDE_M=1414000
# NUM_ORBS=8
# NUM_SATS_PER_ORB=6
# INCLINATION_DEGREE=52

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

ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EARTH_RADIUS = 6378135.0
T = pow((ALTITUDE_M + EARTH_RADIUS) / 21613.546, 3 / 2)
MEAN_MOTION_REV_PER_DAY = round(24 * 60 * 60 / T, 2)

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
distance_file = os.path.join(gen_data, "distance.txt")

sat_info = satgen.read_tles(tle_file)
satellites = sat_info['satellites']
epoch = sat_info['epoch']

isl_list = satgen.read_isls(isl_file, len(satellites))

with open(distance_file, 'w') as f:
    for (a, b) in isl_list:
        sat_distance_m = satgen.distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(epoch))
        f.write(str(sat_distance_m))
        f.write('\n')


