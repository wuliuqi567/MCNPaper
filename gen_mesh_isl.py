import os
import time
import math
import ephem
from astropy.time import Time
from astropy import units as u
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random

from sklearn import preprocessing
# WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html
EARTH_RADIUS = 6378135.0

BASE_NAME = "starlink_550"
NICE_NAME = "Starlink-550"

# STARLINK 550

ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True

################################################################
# The below constants are taken from Starlink's FCC filing as below:
# [1]: https://fcc.report/IBFS/SAT-MOD-20190830-00087
################################################################

MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
ALTITUDE_M = 550000  # Altitude ~550 km

# From https://fcc.report/IBFS/SAT-MOD-20181108-00083/1569860.pdf (minimum angle of elevation: 25 deg)
SATELLITE_CONE_RADIUS_M = 940700

MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))

# ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))

NUM_ORBS = 72
NUM_SATS_PER_ORB = 22
INCLINATION_DEGREE = 53

flows = [1000, 400, 200]  # flow_qos2, flow_qos1 flow_qos0

qos = [[0.33, 0.33, 0.33],
       [0.1, 0.8, 0.1],
       [0.8, 0.1, 0.1]]

gs_stations_file = './gen_data2/starlink_info/gs_with_time.txt'
def get_sats_over_gs_for_time (gs_stations_file):

    gen_data = './gen_data2/'
    sats_info = read_tles(gen_data)
    satellites = sats_info['satellites']
    n_orbits = sats_info['n_orbits']
    n_sats_per_orbit = sats_info['n_sats_per_orbit']
    epoch = sats_info['epoch']

    ground_stations = read_ground_stations_extended(gen_data)


    file_gs = open(gs_stations_file, 'w')

    for i in range(15):
        time = epoch + i * u.minute
        print(time)
        ground_station_satellites_in_range = []
        for ground_station in ground_stations:
            # Find satellites in range
            satellites_in_range = []
            for sid in range(len(satellites)):
                distance_m = distance_m_ground_station_to_satellite(
                    ground_station,
                    satellites[sid],
                    str(epoch),
                    str(time)
                )
                if distance_m <= MAX_GSL_LENGTH_M:
                    satellites_in_range.append((distance_m, sid))
                    # file_gs.write(str(sid))
                    # file_gs.write(',')


            satellites_in_range = sorted(satellites_in_range)
            ground_station_satellites_in_range.append(satellites_in_range)
        # file_gs.write('\n')
    #     ground_station_satellites_in_range_with_time.append(ground_station_satellites_in_range)
    # file_gs.close()

        postions_sat_parsed = get_sat_row_and_column(ground_station_satellites_in_range, sats_info)

        #####################################################################################################
        # sat_two_final_selected = select_minimal_dist(postions_sat_parsed)
        sat_two_final_selected = select_minimal_hop(postions_sat_parsed)

        # sat_two_final_selected = [(num0_plane, num0_sat, dist), (num1_plane, num1_sat, dist)]
        sat_two_final_selected.append((sat_two_final_selected[0][0] * n_sats_per_orbit + sat_two_final_selected[0][1]))
        sat_two_final_selected.append((sat_two_final_selected[1][0] * n_sats_per_orbit + sat_two_final_selected[1][1]))

        sat_connect2gs0_id = sat_two_final_selected[0]
        sat_connect2gs1_id = sat_two_final_selected[1]

        file_gs.write(str(sat_two_final_selected[2]))
        file_gs.write(',')
        file_gs.write(str(sat_two_final_selected[3]))
        file_gs.write('\n')
    file_gs.close()






def get_sats_mesh_net(gen_data, filename_ground_stations_basic_in, filename_ground_stations_out):
    """
    :param gen_data:
    :param filename_ground_stations_basic_in:
    :param filename_ground_stations_out:
    :return:     mesh_info = {
                 "sats_mesh_net_graph": mesh_net,
                 "delay_matrix": delay_matrix,
                 "source2dest_sats": sat_two_final_selected,
                 "mesh_net_pos": pos
                 }
    """

    extend_ground_stations(filename_ground_stations_basic_in, filename_ground_stations_out)
    return gs_selected(gen_data)

def gs_selected(gen_data:str):
    """
    :param gen_data: sats and ground stations' directory
    :return:     mesh_info = {
                 "sats_mesh_net_graph": mesh_net,
                 "delay_matrix": delay_matrix,
                 "source2dest_sats": sat_two_final_selected,
                 "mesh_net_pos": pos
                 }
    """
    #get gs info
    # Dictionary:{
    # "n_orbits": n_orbits,
    # "n_sats_per_orbit": n_sats_per_orbit,
    # "num_of_all_satellite": n_orbits * n_sats_per_orbit,
    # "epoch": epoch,
    # "satellites":satellites
    # }
    sats_info = read_tles(gen_data)
    satellites = sats_info['satellites']
    n_orbits = sats_info['n_orbits']
    n_sats_per_orbit = sats_info['n_sats_per_orbit']
    epoch = sats_info['epoch']
    # time = epoch + 800 * u.s
    ground_stations = read_ground_stations_extended(gen_data)

    sat_data = os.listdir(gen_data)
    file_gs = open(gen_data + sat_data[0] + '/ground_station_satellites_in_range.txt', 'w')

    ground_station_satellites_in_range = []
    for ground_station in ground_stations:
        # Find satellites in range
        satellites_in_range = []
        for sid in range(len(satellites)):
            distance_m = distance_m_ground_station_to_satellite(
                ground_station,
                satellites[sid],
                str(epoch),
                str(epoch)
            )
            if distance_m <= MAX_GSL_LENGTH_M:
                satellites_in_range.append((distance_m, sid))
                file_gs.write(str(sid))
                file_gs.write(',')
                # graph info
                #graphs_sat_net_graph_all_with_only_gsls.add_edge(sid, len(satellites) + ground_station["gid"], weight=distance_m)

        satellites_in_range = sorted(satellites_in_range)
        ground_station_satellites_in_range.append(satellites_in_range)
    file_gs.close()
    # print('sate pos')
    # gs2sat_type = sat_selected(ground_station_satellites_in_range, sat_info)
    postions_sat_parsed = get_sat_row_and_column(ground_station_satellites_in_range, sats_info)

#####################################################################################################
    # sat_two_final_selected = select_minimal_dist(postions_sat_parsed)
    sat_two_final_selected = select_minimal_hop(postions_sat_parsed)
    # sat_two_final_selected ： [(num0_plane, num0_sat, dist), (num1_plane, num1_sat, dist)]
    sat_two_final_selected.append((sat_two_final_selected[0][0] * n_sats_per_orbit + sat_two_final_selected[0][1]))
    sat_two_final_selected.append((sat_two_final_selected[1][0] * n_sats_per_orbit + sat_two_final_selected[1][1]))

    sat_src = sat_two_final_selected[0][0] * n_sats_per_orbit + sat_two_final_selected[0][1]
    sat_des = sat_two_final_selected[1][0] * n_sats_per_orbit + sat_two_final_selected[1][1]

    net_info = get_mesh_net(sat_src, sat_des)

    # add isl according to isl.txt
    mesh_net = nx.Graph()
    mesh_net.add_nodes_from(net_info['sats_in_mesh_list'])
    isl_list = read_isls(gen_data, sats_info['num_of_all_satellite'])
    for (a, b) in isl_list:
        if a in net_info['sats_in_mesh_list'] and b in net_info['sats_in_mesh_list']:
            dist = get_distance_betwenn_adj_sats(satellites[a], satellites[b], sats_info)
            delay = round((dist / 3e8) * 1000, 2)
            mesh_net.add_edge(a, b, weight=delay)

    pos = nx.random_layout(mesh_net)
    nodes_id = list(mesh_net.nodes)

    node_id = 0
    for x in range(len(net_info['orbits_list'])):
        for y in range(len(net_info['n_sat_plane_list'])):
            pos[nodes_id[node_id]][0] = x
            pos[nodes_id[node_id]][1] = -y
            node_id += 1

    #add distance between adjacent-sats

    mesh_info = {
                 "sats_mesh_net_graph": mesh_net,
                 "source2dest_sats": sat_two_final_selected,
                 "mesh_net_pos": pos
                 }
    return mesh_info

def verifiy_routing(isl_list, sat_src, sat_des, sats_info):

    satellites = sats_info['satellites']
    net_info = get_mesh_net(sat_src, sat_des)

    # add isl according to isl.txt
    mesh_net = nx.Graph()
    mesh_net.add_nodes_from(net_info['sats_in_mesh_list'])
    # isl_list = read_isls(gen_data, sats_info['num_of_all_satellite'])
    for (a, b) in isl_list:
        if a in net_info['sats_in_mesh_list'] and b in net_info['sats_in_mesh_list']:
            dist = get_distance_betwenn_adj_sats(satellites[a], satellites[b], sats_info)
            delay = round((dist / 3e8) * 1000, 2)
            mesh_net.add_edge(a, b, weight=delay)

    pos = nx.random_layout(mesh_net)
    nodes_id = list(mesh_net.nodes)

    node_id = 0
    for x in range(len(net_info['orbits_list'])):
        for y in range(len(net_info['n_sat_plane_list'])):
            pos[nodes_id[node_id]][0] = x
            pos[nodes_id[node_id]][1] = -y
            node_id += 1

    minWPath_vs_vt = nx.dijkstra_path(mesh_net, source=sat_src, target=sat_des)
    minWPath_vs_vt_len = nx.dijkstra_path_length(mesh_net, source=sat_src, target=sat_des)

    nx.draw(mesh_net, pos, with_labels=True, font_size=10)
    nx.draw_networkx_nodes(mesh_net, pos, node_color='pink', nodelist=minWPath_vs_vt)

    nx.draw_networkx_nodes(mesh_net, pos, node_color='red', nodelist=[sat_src])
    nx.draw_networkx_nodes(mesh_net, pos, node_color='green', nodelist=[sat_des])

    # node_labels = nx.get_node_attributes(mesh_net, name='pos')
    # nx.draw_networkx_labels(mesh_net, pos, font_size=10, labels=node_labels)
    # edge_labels = nx.get_edge_attributes(mesh_net, 'bandwidth')
    # nx.draw_networkx_edge_labels(mesh_net, pos, font_size=8, edge_labels=edge_labels)

    plt.show()


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



def get_distance_betwenn_adj_sats(sat0_pos, sat1_pos, sats_info):
    """

    :param sat0_pos:
    :param sat1_pos:
    :param sats_info:
    :return: distance_m_between_sats
    """

    satellites = sats_info['satellites']
    n_sats_per_orbit = sats_info['n_sats_per_orbit']

    epoch = sats_info['epoch']
    time = epoch + 0 * u.day
    # parser_sat0_pos = sat0_pos[0] * n_sats_per_orbit + sat0_pos[1]
    # parser_sat1_pos = sat1_pos[0] * n_sats_per_orbit + sat1_pos[1]

    # sat0 = satellites[parser_sat0_pos]
    # sat1 = satellites[parser_sat1_pos]


    return distance_m_between_satellites(sat0_pos, sat1_pos, str(epoch), str(time))

def get_sat_row_and_column(ground_station_satellites_in_range, sats_info):
    """
    :param ground_station_satellites_in_range:  [[], []]
    :param sats_info:
     Dictionary:{
    # "n_orbits": n_orbits,
    # "n_sats_per_orbit": n_sats_per_orbit,
    # "num_of_all_satellite": n_orbits * n_sats_per_orbit,
    # "epoch": epoch,
    # "satellites":satellites
    # }
    :return: #postion_parse=[[(a, b, dist),(c, d, dist)],[(a, b, dist),(c, d, dist)]]
    """
    postions_sat_parsed = []

    n_sats_per_orbit = sats_info['n_sats_per_orbit']

    for ground_station in ground_station_satellites_in_range:
        # parse num_plan and num_sat_in_plane
        # mini_dist = ground_station[0][0]
        temp = []
        for dist, id in ground_station:
            # if dist <= 1.3 * mini_dist:
            # if dist <= mini_dist:
            temp.append((id // n_sats_per_orbit, id % n_sats_per_orbit, dist))
        postions_sat_parsed.append(temp)
    # temp = []
    # id = ground_station_satellites_in_range[0][0][1]
    # dist = ground_station_satellites_in_range[0][0][0]
    # temp.append((id // n_sats_per_orbit, id % n_sats_per_orbit, dist))
    # postions_sat_parsed.append(temp)
    #
    # temp = []
    # id = ground_station_satellites_in_range[1][0][1]
    # dist = ground_station_satellites_in_range[1][0][0]
    # temp.append((id // n_sats_per_orbit, id % n_sats_per_orbit, dist))

    # postions_sat_parsed.append(temp)


    return postions_sat_parsed

def select_minimal_hop(postions_sat_parsed):
    sat_two_final_selected = []

    minimal_mesh_grid = 65535
    for a, b, dist0 in postions_sat_parsed[0]:
        for c, d, dist1 in postions_sat_parsed[1]:
            # if a == c:
            #     rows_nums = 1
            # else:
            #     rows_nums = get_positive_int((a - c))
            rows_nums = min(get_positive_int((a - c)), 72-get_positive_int((a - c)))
            # if b == d:
            #     column = 1
            # else:
            #     column = get_positive_int((b - d))
            column = min(get_positive_int((b - d)), 22-get_positive_int((b - d)))
            num_grid = rows_nums + column
            if get_positive_int(num_grid) < minimal_mesh_grid:
                sat0_orbit = a
                sat0_num_of_plane = b
                final_row = rows_nums
                sat0_to_gs_dist = dist0

                sat1_orbit = c
                sat1_num_of_plane = d
                sat1_to_gs_dist = dist1
                minimal_mesh_grid = get_positive_int(num_grid)
            elif num_grid == minimal_mesh_grid:
                if rows_nums > final_row:
                #if dist0 + dist1 < sat0_to_gs_dist + sat1_to_gs_dist:
                    sat0_orbit = a
                    sat0_num_of_plane = b
                    final_row = rows_nums
                    sat0_to_gs_dist = dist0

                    sat1_orbit = c
                    sat1_num_of_plane = d
                    sat1_to_gs_dist = dist1

    sat_two_final_selected.append((sat0_orbit, sat0_num_of_plane, sat0_to_gs_dist))
    sat_two_final_selected.append((sat1_orbit, sat1_num_of_plane, sat1_to_gs_dist))

    return sat_two_final_selected

def select_minimal_dist(postions_sat_parsed):
    sat_two_final_selected = []

    sat_two_final_selected.append(postions_sat_parsed[0][1])
    sat_two_final_selected.append(postions_sat_parsed[1][1])
    return sat_two_final_selected

def get_positive_int(ne_or_po_num:int):
    if ne_or_po_num < 0:
        ne_or_po_num = - ne_or_po_num
    return ne_or_po_num

def sat_selected(ground_station_satellites_in_range:list, sat_info:dict):
    """

    :param ground_station_satellites_in_range:
    :param sat_info:
    :return: gs2sat_type = {"gs1":{"up":[(distance, sid),(distance, sid)], "down":[(distance, sid),(distance, sid)]}, "gs1":{"up":[(distance, sid),(distance, sid)]}
    """


    satellites = sat_info['satellites']
    epoch = sat_info['epoch']
    time = epoch + 0 * u.day
    d = ephem.Date(str(time))

    gs2sat_type = {}
    # {"gs1":{"up":[(distance, sid),(distance, sid)], "down":[(distance, sid),(distance, sid)]}, "gs1":{"up":[(distance, sid),(distance, sid)]},
    gs2sat_type.setdefault(0, {'up':[], 'down':[]})
    gs2sat_type.setdefault(1, {'up':[], 'down':[]})
    # print(gs2sat_type)

    for id, ground_station in enumerate(ground_station_satellites_in_range):
        for distance, sat_id in ground_station:
            #lable sat_id is ascending or descending

            satellite = satellites[sat_id]
            satellite.compute(d)

            if satellite.ra >= ephem.degrees(str(satellite.raan-90)) and satellite.ra < ephem.degrees(str(satellite.raan-90)):
                #
                print('ascending', satellite.ra)
                gs2sat_type[id]['up'].append((distance, sat_id))
            else:
                print('descending', satellite.ra)
                #
                gs2sat_type[id]['down'].append((distance, sat_id))
        # print('-------------------')
        # print(gs2sat_type)
        # print('------------')
        for id in gs2sat_type:
            for id2 in gs2sat_type[id]:
                sorted(gs2sat_type[id][id2])
                # print(id2)
                # print(gs2sat_type[id][id2])
        """
        # print('%s %s' % (satellite.sublong, satellite.sublat))  # 卫星的星下点
                    # print(float(satellite.ra))  # 赤经
                    # print(repr(satellite.ra))

                    # print(satellite.ra / pi * 180.0)
                    # print(satellite.ra)   # 赤经
                    # print(satellite.dec)  # 赤纬
                    # # print(satellite.elevation)
                    # print(satellite.raan)
                    # print(satellite.inc)
                    # print(satellite.n)
                    #
                # EarthSatellite elements of man-made satellites:
                #
                # epoch — Reference epoch
                # n — Mean motion, in revolutions per day
                # inc — Inclination (°)
                # raan — Right Ascension of ascending node (°)
                # e — Eccentricity
                # ap — Argument of perigee at epoch (°)
                # M — Mean anomaly from perigee at epoch (°)
                # decay — Orbit decay rate in revolutions per day, per day
                # drag — Object drag coefficient in per earth radii
                # orbit — Integer orbit number of epoch
                    # print('------------')
        """
        return  gs2sat_type

def read_tles(gen_data):
    """
    Read a constellation of satellites from the TLES file.

    :param filename_tles:                    Filename of the TLES (typically /path/to/tles.txt)

    :return: Dictionary:{
        "n_orbits": n_orbits,
        "n_sats_per_orbit": n_sats_per_orbit,
        "num_of_all_satellite": n_orbits * n_sats_per_orbit,
        "epoch": epoch,
        "satellites":satellites
        }
    """

    sat_data = os.listdir(gen_data)
    satellites = []
    with open(gen_data + sat_data[0] +'/tles.txt', 'r') as f:
        n_orbits, n_sats_per_orbit = [int(n) for n in f.readline().split()]
        universal_epoch = None

        i = 0
        for tles_line_1 in f:
            tles_line_2 = f.readline()
            tles_line_3 = f.readline()

            # Retrieve name and identifier
            name = tles_line_1
            sid = int(name.split()[1])
            if sid != i:
                raise ValueError("Satellite identifier is not increasing by one each line")
            i += 1

            # Fetch and check the epoch from the TLES data
            # In the TLE, the epoch is given with a Julian data of yyddd.fraction
            # ddd is actually one-based, meaning e.g. 18001 is 1st of January, or 2018-01-01 00:00.
            # As such, to convert it to Astropy Time, we add (ddd - 1) days to it
            # See also: https://www.celestrak.com/columns/v04n03/#FAQ04
            epoch_year = tles_line_2[18:20]
            epoch_day = float(tles_line_2[20:32])
            epoch = Time("20" + epoch_year + "-01-01 00:00:00", scale="tdb") + (epoch_day - 1) * u.day
            if universal_epoch is None:
                universal_epoch = epoch
            if epoch != universal_epoch:
                raise ValueError("The epoch of all TLES must be the same")

            # Finally, store the satellite information
            satellites.append(ephem.readtle(tles_line_1, tles_line_2, tles_line_3))

    return {
        "n_orbits": n_orbits,
        "n_sats_per_orbit": n_sats_per_orbit,
        "num_of_all_satellite": n_orbits * n_sats_per_orbit,
        "epoch": epoch,
        "satellites":satellites
    }

def read_ground_stations_extended(gen_filename):
    """
    Reads ground stations from the input file.

    :param filename_ground_stations_basic: Filename of ground stations basic (typically /path/to/ground_stations.txt)

    :return: List of ground stations
    """

    constellation_info_list = os.listdir(gen_filename)

    ground_stations_extended = []
    gid = 0

    for i in constellation_info_list:
        with open(gen_filename + i +'/ground_stations.txt', 'r') as f:
            for line in f:
                split = line.split(',')
                if len(split) != 8:
                    raise ValueError("Extended ground station file has 8 columns: " + line)
                if int(split[0]) != gid:
                    raise ValueError("Ground station id must increment each line")
                ground_station_basic = {
                    "gid": gid,
                    "name": split[1],
                    "latitude_degrees_str": split[2],
                    "longitude_degrees_str": split[3],
                    "elevation_m_float": float(split[4]),
                    "cartesian_x": float(split[5]),
                    "cartesian_y": float(split[6]),
                    "cartesian_z": float(split[7]),
                }
                ground_stations_extended.append(ground_station_basic)
                gid += 1
    return ground_stations_extended

def read_isls(gen_filename, num_satellites):
    """
    Read ISLs file into a list of undirected edges

    :param filename_isls:  Filename of ISLs (typically /path/to/isls.txt)
    :param num_satellites: Number of satellites (to verify indices)

    :return: List of all undirected ISL edges
    """
    isls_list = []

    info = os.listdir(gen_filename)

    with open(gen_filename + info[0] + '/isls.txt', 'r') as f:
        isls_set = set()
        for line in f:
            line_spl = line.split()
            a = int(line_spl[0])
            b = int(line_spl[1])

            # Verify the input
            if a >= num_satellites:
                raise ValueError("Satellite does not exist: %d" % a)
            if b >= num_satellites:
                raise ValueError("Satellite does not exist: %d" % b)
            if b <= a:
                raise ValueError("The second satellite index must be strictly larger than the first")
            if (a, b) in isls_set:
                raise ValueError("Duplicate ISL: (%d, %d)" % (a, b))
            isls_set.add((a, b))

            # Finally add it to the list
            isls_list.append((a, b))

    return isls_list

def distance_m_between_satellites(sat1, sat2, epoch_str, date_str):
    """
    Computes the straight distance between two satellites in meters.

    :param sat1:       The first satellite
    :param sat2:       The other satellite
    :param epoch_str:  Epoch time of the observer (string)
    :param date_str:   The time instant when the distance should be measured (string)

    :return: The distance between the satellites in meters
    """

    # Create an observer somewhere on the planet
    observer = ephem.Observer()
    observer.epoch = epoch_str
    observer.date = date_str
    observer.lat = 0
    observer.lon = 0
    observer.elevation = 0

    # Calculate the relative location of the satellites to this observer
    sat1.compute(observer)
    sat2.compute(observer)

    # Calculate the angle observed by the observer to the satellites (this is done because the .compute() calls earlier)
    angle_radians = float(repr(ephem.separation(sat1, sat2)))

    # Now we have a triangle with three knows:
    # (1) a = sat1.range (distance observer to satellite 1)
    # (2) b = sat2.range (distance observer to satellite 2)
    # (3) C = angle (the angle at the observer point within the triangle)
    #
    # Using the formula:
    # c^2 = a^2 + b^2 - 2 * a * b * cos(C)
    #
    # This gives us side c, the distance between the two satellites
    return math.sqrt(sat1.range ** 2 + sat2.range ** 2 - (2 * sat1.range * sat2.range * math.cos(angle_radians)))

def distance_m_ground_station_to_satellite(ground_station, satellite, epoch_str, date_str):
    """
    Computes the straight distance between a ground station and a satellite in meters

    :param ground_station:  The ground station
    :param satellite:       The satellite
    :param epoch_str:       Epoch time of the observer (ground station) (string)
    :param date_str:        The time instant when the distance should be measured (string)

    :return: The distance between the ground station and the satellite in meters
    """

    # Create an observer on the planet where the ground station is
    observer = ephem.Observer()
    observer.epoch = epoch_str
    observer.date = date_str
    observer.lat = str(ground_station["latitude_degrees_str"])   # Very important: string argument is in degrees.
    observer.lon = str(ground_station["longitude_degrees_str"])  # DO NOT pass a float as it is interpreted as radians
    observer.elevation = ground_station["elevation_m_float"]

    # Compute distance from satellite to observer
    satellite.compute(observer)

    # Return distance
    return satellite.range

def read_ground_stations_basic(filename_ground_stations_basic):
    """
    Reads ground stations from the input file.

    :param filename_ground_stations_basic: Filename of ground stations basic (typically /path/to/ground_stations.txt)

    :return: List of ground stations
    """
    ground_stations_basic = []
    gid = 0
    with open(filename_ground_stations_basic, 'r') as f:
        for line in f:
            split = line.split(',')
            if len(split) != 5:
                raise ValueError("Basic ground station file has 5 columns")
            if int(split[0]) != gid:
                raise ValueError("Ground station id must increment each line")
            ground_station_basic = {
                "gid": gid,
                "name": split[1],
                "latitude_degrees_str": split[2],
                "longitude_degrees_str": split[3],
                "elevation_m_float": float(split[4]),
            }
            ground_stations_basic.append(ground_station_basic)
            gid += 1
    return ground_stations_basic

def extend_ground_stations(filename_ground_stations_basic_in, filename_ground_stations_out):
    ground_stations = read_ground_stations_basic(filename_ground_stations_basic_in)
    with open(filename_ground_stations_out, "w+") as f_out:
        for ground_station in ground_stations:
            cartesian = geodetic2cartesian(
                float(ground_station["latitude_degrees_str"]),
                float(ground_station["longitude_degrees_str"]),
                ground_station["elevation_m_float"]
            )
            f_out.write(
                "%d,%s,%f,%f,%f,%f,%f,%f\n" % (
                    ground_station["gid"],
                    ground_station["name"],
                    float(ground_station["latitude_degrees_str"]),
                    float(ground_station["longitude_degrees_str"]),
                    ground_station["elevation_m_float"],
                    cartesian[0],
                    cartesian[1],
                    cartesian[2]
                )
            )

def geodetic2cartesian(lat_degrees, lon_degrees, ele_m):
    """
    Compute geodetic coordinates (latitude, longitude, elevation) to Cartesian coordinates.

    :param lat_degrees: Latitude in degrees (float)
    :param lon_degrees: Longitude in degrees (float)
    :param ele_m:  Elevation in meters

    :return: Cartesian coordinate as 3-tuple of (x, y, z)
    """

    #
    # Adapted from: https://github.com/andykee/pygeodesy/blob/master/pygeodesy/transform.py
    #

    # WGS72 value,
    # Source: https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html
    a = 6378135.0

    # Ellipsoid flattening factor; WGS72 value
    # Taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html
    f = 1.0 / 298.26

    # First numerical eccentricity of ellipsoid
    e = math.sqrt(2.0 * f - f * f)
    lat = lat_degrees * (math.pi / 180.0)
    lon = lon_degrees * (math.pi / 180.0)

    # Radius of curvature in the prime vertical of the surface of the geodetic ellipsoid
    v = a / math.sqrt(1.0 - e * e * math.sin(lat) * math.sin(lat))

    x = (v + ele_m) * math.cos(lat) * math.cos(lon)
    y = (v + ele_m) * math.cos(lat) * math.sin(lon)
    z = (v * (1.0 - e * e) + ele_m) * math.sin(lat)

    return x, y, z

def get_qoslinks_weight(qosRank, delay_matrix, bandwidth, packet_lr):
    """
    :param delay_matrix:
    :return: qoslink_weith = [qos_list]
    """
    """
    The bandwidth of the different links is the randomly chosen from 500 KBps to 2 MBps.
    The packet loss rates of the different links are randomly chosen from 0.2 to 1.5 percent
    """

    # delay_matrix = maxminnormalization(delay_matrix)
    # bandwidth = maxminnormalization(bandwidth)
    # packet_lr = maxminnormalization(packet_lr)

    # gen adj yinzi
    Qb = round(min(bandwidth) / sum(bandwidth), 4)
    Qd = round(min(delay_matrix) / sum(delay_matrix), 4)
    Ql = round(min(packet_lr) / sum(packet_lr), 6)

    qos_links_weight =[]
    if qosRank == 2:
        qos2link_weith = []
        qos2 = qos[2]
        for index in range(len(bandwidth)):
            weight = round(qos2[0] * Qb * bandwidth[index] - qos2[1] * Qd * delay_matrix[index] - qos2[2] * Ql * packet_lr[index], 2)
            # if weight <= 0:
            #     weight = 0.01
            qos2link_weith.append(weight)
        min_link_weidth = get_positive_int(min(qos2link_weith)) + 0.01

        for idx in range(len(qos2link_weith)):
            qos2link_weith[idx] += min_link_weidth

        return qos2link_weith
        # qos_links_weight.append(maxminnormalization(qos2link_weith))

    if qosRank == 1:
        qos0link_weith = []
        qos0 = qos[0]
        for index in range(len(bandwidth)):
            weight = round(qos0[0] * Qb * bandwidth[index] - qos0[1] * Qd * delay_matrix[index] - qos0[2] * Ql * packet_lr[index], 2)
            # if weight <= 0:
            #     weight = 0.01
            qos0link_weith.append(weight)
        min_link_weidth = get_positive_int(min(qos0link_weith)) + 0.01
        for idx in range(len(qos0link_weith)):
            qos0link_weith[idx] += min_link_weidth
        return qos0link_weith
        # qos_links_weight.append(maxminnormalization(qos0link_weith))

    if qosRank == 0:
        qos1link_weith = []
        qos1 = qos[1]
        for index in range(len(bandwidth)):
            weight = round(qos1[0] * Qb * bandwidth[index] - qos1[1] * Qd * delay_matrix[index] - qos1[2] * Ql * packet_lr[index], 2)
            # if weight <= 0:
            #     weight = 0.01
            qos1link_weith.append(weight)
        min_link_weidth = get_positive_int(min(qos1link_weith)) + 0.01
        for idx in range(len(qos1link_weith)):
            qos1link_weith[idx] += min_link_weidth
        return qos1link_weith
        # qos_links_weight.append(maxminnormalization(qos1link_weith))

    # return qos_links_weight

def maxminnormalization(x:list):
    Max = max(x)
    Min = min(x)
    new_x = []
    for i in range(len(x)):
        new_x.append(round((x[i] - Min) / (Max - Min), 2))
        if new_x[i] == 0:
            new_x[i] = 0.01
        # x[i] = round(1/x[i], 2)
    return new_x


if __name__=='__main__':

    gen_data = './gen_data2/'
    filename_ground_stations_basic_in = './gen_data2/starlink_info/ground_stations_basic.txt'
    filename_ground_stations_out = './gen_data2/starlink_info/ground_stations.txt'

    mesh_info = get_sats_mesh_net(gen_data, filename_ground_stations_basic_in, filename_ground_stations_out)
    """
        mesh_info = {
                 "sats_mesh_net_graph": mesh_net,
                 "delay_matrix": delay_matrix,
                 "source2dest_sats": sat_two_final_selected,
                 "mesh_net_pos": pos
                 }
    """

    mesh_net = mesh_info['sats_mesh_net_graph']
    source2dest_sats = mesh_info["source2dest_sats"]
    pos = mesh_info["mesh_net_pos"]

    nx.draw(mesh_net, pos)
    # nx.draw_networkx_nodes(mesh_net, pos, node_color=node_color[i], nodelist=minWPath_vs_vt)

    nx.draw_networkx_nodes(mesh_net, pos, node_color='red', nodelist=[source2dest_sats[2]])
    nx.draw_networkx_nodes(mesh_net, pos, node_color='yellow', nodelist=[source2dest_sats[3]])

    node_labels = nx.get_node_attributes(mesh_net, name='pos')
    nx.draw_networkx_labels(mesh_net, pos, font_size=10, labels=node_labels)
    edge_labels = nx.get_edge_attributes(mesh_net, 'delay')
    nx.draw_networkx_edge_labels(mesh_net, pos, font_size=8, edge_labels=edge_labels)

    plt.show()
