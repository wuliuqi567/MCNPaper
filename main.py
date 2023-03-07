
import os
from gen_satellites_info import gen_info


# read configure file
configure_dict = {}
# {'satellite_network_dir_and_name': 'starlink',
# 'ALTITUDE_M': '550000',
# 'Elevation_angle': '20',
# 'NUM_ORBS': '72',
# 'NUM_SATS_PER_ORB': '22',
# 'INCLINATION_DEGREE': '53',
# 'simulation_end_time_s': '600',
# 'Time_step_ms': '100',
# 'isl_data_rate_megabit_per_s': '10.0',
# 'gsl_data_rate_megabit_per_s': '10.0',
# 'isl_max_queue_size_pkts': '100',
# 'gsl_max_queue_size_pkts': '100',
# }
with open('config_mininet.csv', 'r') as config_file:
    lines = config_file.readlines()
    for line in lines:
        line = line.strip().split('=')
        configure_dict.setdefault(line[0], line[1])
    # print(configure_dict)xc

if not os.path.exists(configure_dict['satellite_network_dir']):
    #gen_info(output_generated_data_dir, simulation_end_time_s, time_step_ms,
    #             BASE_NAME, ALTITUDE_M, Elevation_angle, NUM_ORBS, NUM_SATS_PER_ORB, INCLINATION_DEGREE):
    gen_info(configure_dict['satellite_network_dir_and_name'], int(configure_dict['simulation_end_time_s']),
             int(configure_dict['Time_step_ms']), configure_dict['satellite_network_dir_and_name'],
             int(configure_dict['ALTITUDE_M']), int(configure_dict['Elevation_angle']), int(configure_dict['NUM_ORBS']),
             int(configure_dict['NUM_SATS_PER_ORB']), int(configure_dict['INCLINATION_DEGREE']))

else:
    print('baibai')

