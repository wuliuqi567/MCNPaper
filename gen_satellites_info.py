import sys
import math
from main_helper import MainHelper

EARTH_RADIUS = 6378135.0

ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True

################################################################
# The below constants are taken from Starlink's FCC filing as below:
# [1]: https://fcc.report/IBFS/SAT-MOD-20190830-00087
################################################################
# BASE_NAME = "starlink_550"
# ALTITUDE_M = 550000  # Altitude ~550 km
# Elevation_angle = 25
# NUM_ORBS = 72
# NUM_SATS_PER_ORB = 22
# INCLINATION_DEGREE = 53


def gen_info(output_generated_data_dir, simulation_end_time_s, time_step_ms,
             BASE_NAME, ALTITUDE_M, Elevation_angle, NUM_ORBS, NUM_SATS_PER_ORB, INCLINATION_DEGREE):

    # output_generated_data_dir = '../starlink/'  # Final directory in which the result will be placed
    # simulation_end_time_s = 1500
    # time_step_ms = 10000
    gs_selection = "twostation" # ground_stations_{top_100, paris_moscow_grid}

    isl_selection = 'isls_plus_grid'  # isls_{none, plus_grid}
    dynamic_state_algorithm = 'algorithm_free_one_only_over_isls' # algorithm_{free_one_only_{gs_relays,_over_isls}, paired_many_only_over_isls}
    num_threads = 1

    T = pow((ALTITUDE_M + EARTH_RADIUS) / 21613.546, 3 / 2)
    MEAN_MOTION_REV_PER_DAY = round(24 * 60 * 60 / T, 2)

    SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(Elevation_angle))
    MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))

    # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
    MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))

    main_helper = MainHelper(
        BASE_NAME,
        BASE_NAME,
        ECCENTRICITY,
        ARG_OF_PERIGEE_DEGREE,
        PHASE_DIFF,
        MEAN_MOTION_REV_PER_DAY,
        ALTITUDE_M,
        MAX_GSL_LENGTH_M,
        MAX_ISL_LENGTH_M,
        NUM_ORBS,
        NUM_SATS_PER_ORB,
        INCLINATION_DEGREE,
    )

    main_helper.calculate(
        output_generated_data_dir,
        simulation_end_time_s,
        time_step_ms,
        isl_selection,
        gs_selection,
        dynamic_state_algorithm,
        num_threads,
    )
#python main_starlink_550.py 3600 600000 isls_plus_grid ground_stations_top_100 algorithm_free_one_only_over_isls 30
