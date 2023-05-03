import os
import math
import ephem
from astropy.time import Time
from astropy import units as u
import networkx as nx
import time
from gen_satellites_info import gen_info
from gen_mesh_isl import verifiy_routing

EARTH_RADIUS = 6378135.0
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
light_speed = 3e8

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
    # print(configure_dict)

if not os.path.exists(configure_dict['satellite_network_dir_and_name']):
    gen_info(configure_dict['satellite_network_dir_and_name'], int(configure_dict['simulation_end_time_s']),
             int(configure_dict['Time_step_ms']), configure_dict['satellite_network_dir_and_name'],
             int(configure_dict['ALTITUDE_M']), int(configure_dict['Elevation_angle']), int(configure_dict['NUM_ORBS']),
             int(configure_dict['NUM_SATS_PER_ORB']), int(configure_dict['INCLINATION_DEGREE']))
else:
    print('constellation {} has already'.format(configure_dict['satellite_network_dir_and_name']))

ALTITUDE_M = int(configure_dict['ALTITUDE_M'])
NUM_ORBS = int(configure_dict['NUM_ORBS'])
NUM_SATS_PER_ORB = int(configure_dict['NUM_SATS_PER_ORB'])
elevation_angle = int(configure_dict['Elevation_angle'])
simulation_end_time_s = int(configure_dict['simulation_end_time_s'])
time_step_ms = int(configure_dict['Time_step_ms'])
gen_data = configure_dict['satellite_network_dir_and_name']

isl_data_rate_megabit_per_s = float(configure_dict['isl_data_rate_megabit_per_s'])
gsl_data_rate_megabit_per_s = float(configure_dict['gsl_data_rate_megabit_per_s'])
isl_max_queue_size_pkts = int(configure_dict['isl_max_queue_size_pkts'])
gsl_max_queue_size_pkts = int(configure_dict['gsl_max_queue_size_pkts'])

satellite_network_dir = configure_dict['satellite_network_dir_and_name']
satellite_network_routes_dir = os.path.join(satellite_network_dir, satellite_network_dir, ) \
                               + '_isls_plus_grid_twostation_algorithm_free_one_only_over_isls/dynamic_state_{}ms_for_{}s'.format(
    time_step_ms, simulation_end_time_s)

# Considering an elevation angle of 30 degrees; possible values [1]: 20(min)/30/35/45
SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(elevation_angle))

MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))

graphs_sat_net_graph_only_satellites_with_isls = nx.Graph()
graphs_sat_net_graph_all_with_only_gsls = nx.Graph()


enable_verbose_logs = True
# network Graphs

# interfaces info
# self.total_num_isls = 0
# self.num_isls_per_sat = []
# self.sat_neighbor_to_if = {}


def OGDRoutingTest(satellites, ground_station_satellites_in_range, graphs_sat_net):

    access_sat = calculate_two_gs_routing(satellites, ground_station_satellites_in_range, 0, 1)
    minWPath_vs_vt = nx.dijkstra_path(graphs_sat_net, source=access_sat[0], target=access_sat[1])
    return minWPath_vs_vt


# [[(614842.25, 1521, 1), (683955.5, 1085, 0), (687683.5625, 1064, 0), (693692.875, 1543, 1),
def OGDRouting(satellites, ground_station_satellites_in_range):
    # divide satellites into two types
    gs_size = len(ground_station_satellites_in_range)
    route = []
    for gs1_index in range(0, gs_size):
        gs2_index = gs_size - 1
        while gs2_index > gs1_index:
            temp = calculate_two_gs_routing(satellites, ground_station_satellites_in_range, gs1_index, gs2_index)
            temp.insert(0, len(satellites) + gs1_index)
            temp.append(len(satellites) + gs2_index)
            route.append(temp)
            gs2_index -= 1
    return route


def calculate_two_gs_routing(satellites, ground_station_satellites_in_range, src_gs_index, des_gs_index):
    src_gs_sats_in_range = ground_station_satellites_in_range[src_gs_index]
    des_gs_sats_in_range = ground_station_satellites_in_range[des_gs_index]
    min_hop_sum = 65535
    select_src_sat_id = 0
    select_des_sat_id = 0
    # dir_hop_horizontal = 0
    # dir_hop_vertical = 0
    hop_horizontal_last = 0
    sel_src_type = 2  # 0 descending 1 ascending  2 error
    sel_des_type = 2

    for dis_sat_to_src, src_sat_id, src_type in src_gs_sats_in_range:
        for dis_sat_to_des, des_sat_id, des_type in des_gs_sats_in_range:

            n_src = src_sat_id // 22
            n_des = des_sat_id // 22

            m_src = src_sat_id % 22
            m_des = des_sat_id % 22

            hop_horizontal = min(abs(n_src - n_des), 72 - abs(n_src - n_des))
            hop_vertical = min(abs(m_src - m_des), 22 - abs(m_src - m_des))

            hop_sum = hop_horizontal + hop_vertical
            # print('src  des  hop', src_sat_id, des_sat_id, hop_sum)
            if hop_sum < min_hop_sum:
                min_hop_sum = hop_sum
                select_src_sat_id = src_sat_id
                select_des_sat_id = des_sat_id
                hop_horizontal_last = hop_horizontal
                sel_src_type = src_type
                sel_des_type = des_type

            elif hop_sum == min_hop_sum:
                if hop_horizontal > hop_horizontal_last:
                    min_hop_sum = hop_sum
                    select_src_sat_id = src_sat_id
                    select_des_sat_id = des_sat_id
                    hop_horizontal_last = hop_horizontal
                    sel_src_type = src_type
                    sel_des_type = des_type

    return [select_src_sat_id, select_des_sat_id]

    # print('select_src_sat_id  select_des_sat_id hop', select_src_sat_id, select_des_sat_id, min_hop_sum)
    # if sel_src_type == sel_des_type:
    #     # return orbit_gird_routing_sametype(satellites, select_src_sat_id, select_des_sat_id)
    #     return OGDRouting_A2A(satellites, select_src_sat_id, select_des_sat_id)
    # else:
    #     return OGDRouting_A2D(satellites, select_src_sat_id, select_des_sat_id)


def OGDRouting_A2A(satellites, select_src_sat_id, select_des_sat_id):
    route_from_s = [select_src_sat_id]
    route_from_d = [select_des_sat_id]

    n_src = select_src_sat_id // 22
    n_des = select_des_sat_id // 22
    hop_h = n_src - n_des
    if -72 / 2 <= hop_h < 0 or hop_h > 72 / 2:
        dir_hop_horizontal = 1  # right 顺时针减 逆时针加
    elif hop_h < -72 / 2 or 0 < hop_h <= 72 / 2:
        dir_hop_horizontal = -1  # left
    else:
        dir_hop_horizontal = 0  # no movement

    m_src = select_src_sat_id % 22
    m_des = select_des_sat_id % 22
    hop_v = m_src - m_des
    if -22 / 2 <= hop_v < 0 or hop_v > 22 / 2:
        dir_hop_vertical = 1  # up 顺时针加 逆时针减
    elif hop_v < -22 / 2 or 0 < hop_v <= 22 / 2:
        dir_hop_vertical = -1  # down
    else:
        dir_hop_vertical = 0  # no movement

    while len(route_from_s) > 0 and len(route_from_d) > 0 and route_from_s[-1] // 22 != route_from_d[
        -1] // 22:  # 直到 n_src == n_des

        curr_node_src = route_from_s[-1]
        curr_node_des = route_from_d[-1]
        next_node_n_src = Get_N_of_plane(curr_node_src // 22, dir_hop_horizontal, '+')
        next_node_n_des = Get_N_of_plane(curr_node_des // 22, dir_hop_horizontal, '-')
        next_node_src = next_node_n_src * 22 + curr_node_src % 22
        next_node_des = next_node_n_des * 22 + curr_node_des % 22

        reward_s = abs(satellites[curr_node_src].sublat) + abs(satellites[next_node_src].sublat)
        reward_d = abs(satellites[curr_node_des].sublat) + abs(satellites[next_node_des].sublat)

        if reward_s > reward_d:

            if LinkConnected(curr_node_src, next_node_src):  # 向前正常
                route_from_s.append(next_node_src)
            else:  # 向前链路异常
                next_node_m_src = Get_M_of_plane(curr_node_src % 22, dir_hop_vertical, "+")  # 向上查找
                next_node_src = (curr_node_src // 22) * 22 + next_node_m_src
                if LinkConnected(curr_node_src, next_node_src):  # 向上链路正常
                    route_from_s.append(next_node_src)
                else:  # 向上链路异常 回退
                    route_from_s = HorizontalBacktrackingForA2A(route_from_s, '+', dir_hop_vertical)  # 无路可退返回空list
        else:
            if LinkConnected(curr_node_des, next_node_des):
                route_from_d.append(next_node_des)
            else:  # 向后链路异常
                next_node_m_des = Get_M_of_plane(curr_node_des % 22, dir_hop_vertical, '-')
                next_node_des = (curr_node_des // 22) * 22 + next_node_m_des
                if LinkConnected(curr_node_des, next_node_des):  # 向下的链路正常
                    route_from_d.append(next_node_des)
                else:  # 向下的链路异常 回退
                    route_from_d = HorizontalBacktrackingForA2A(route_from_d, '-', dir_hop_vertical)  # 无路可退返回空list

    if len(route_from_s) < 1 or len(route_from_d) < 1:
        return []
    elif route_from_s[-1] // 22 != route_from_d[-1] // 22:
        return []

    nodes = VerticalLinksCheckForA2A(route_from_s, route_from_d, dir_hop_vertical)
    if len(nodes) < 1:  # 链路故障 返回空为链路故障
        nodes, route_from_s, route_from_d = VerticalSearchForA2A(route_from_s, route_from_d, dir_hop_horizontal,
                                                                 dir_hop_vertical)
        if len(route_from_s) > 0:  # 找到了路
            route_from_s += nodes[1:-1]
        else:
            return  # 无路可走
    else:  # 链路正常，路由计算完毕
        route_from_s += nodes[1:-1]

    route_from_d.reverse()
    return route_from_s + route_from_d


def HorizontalBacktrackingForA2A(route, type_from_s_or_d, dir_hop_vertical):
    """
    向前向上链路故障时候进入，用于链路回退/回溯
    解释说明：
    for example:
    route = [1,2,3,4,5]
    temp_node->5
    curr_node->4
    5节点向前向上链路故障，回到4节点进行链路探查
    如果 temp_node 和 curr_node 处于同一轨道内编号，那么从curr_node向上探查
        如果 向上探查到链路正常，将其加入路径
        否则 链路故障 再次回退
    否则 再次回退
    :param route:
    :param type_from_s_or_d: 如果是 route = route_from_s 应该设为 '+'，否则 '-' 表示路由寻找的方向
    :param dir_hop_vertical:
    :return: 无路可退返回空list
    """
    if len(route) == 1:  # 向前向上链路故障才进来，所以如果是从源/端点开始向上先前链路故障即为无路
        return []
    temp_node = route.pop()
    curr_node = route[-1]

    if temp_node % 22 == curr_node % 22:  # 横向退 m = m
        next_node_m_src = Get_M_of_plane(curr_node % 22, dir_hop_vertical, type_from_s_or_d)  # if s '+' else '-'
        next_node_src = (curr_node // 22) * 22 + next_node_m_src
        if LinkConnected(curr_node, next_node_src):  # 向上检查
            route.append(next_node_src)
        else:  # 链路不存在 再次回退
            route = HorizontalBacktrackingForA2A(route, type_from_s_or_d, dir_hop_vertical)
    else:  # 纵向退
        route = HorizontalBacktrackingForA2A(route, type_from_s_or_d, dir_hop_vertical)

    return route


def VerticalSearchForA2A(route_from_s, route_from_d, dir_hop_horizontal, dir_hop_vertical):
    """
    在左右范围内上下寻找，
    - - - - - - - - - 6 d
    - - - - - - 2 3 4 1
    - - - - - - v -- -
    - - - - 1 3 5 - - -
    s 2 4 5 4
    即，在1 与 1 之间找路 采用向左向右循环找，找到一个即可
    :param route:
    :param route_from_s:
    :param route_from_d:
    :param dir_hop_horizontal:
    :param dir_hop_vertical:
    :return: new_route_from_s = []代表无路，
    """

    Rs_src = route_from_s
    Rd_src = route_from_d
    Rs_des = route_from_s
    Rd_des = route_from_d

    new_route_from_s = []
    new_route_from_d = []

    nodes = []
    flag_1 = True
    flag_2 = True
    while flag_1 or flag_2:

        next_n_src = Get_N_of_plane(Rs_src[-1] // 22, dir_hop_horizontal, '+')
        next_node_src = next_n_src * 22 + Rs_src[-1] % 22

        if flag_1 and len(Rd_src) > 1 and Rd_src[-1] % 22 == Rd_src[-2] % 22 and LinkConnected(Rs_src[-1],
                                                                                               next_node_src):  # 向前链路正常
            Rd_src.pop()
            Rs_src.append(next_node_src)
            nodes = VerticalLinksCheckForA2A(Rs_src, Rd_src, dir_hop_vertical)
            if len(nodes) > 0:  # 垂直链路正常
                new_route_from_s = Rs_src
                new_route_from_d = Rd_src
                break
        else:
            flag_1 = False  # 向前链路异常 向前结束

        next_n_des = Get_N_of_plane(Rd_des[-1] // 22, dir_hop_horizontal, '-')
        next_node_des = next_n_des * 22 + Rd_des[-1] % 22

        if flag_2 and len(Rs_des) > 1 and Rs_des[-1] % 22 == Rs_des[-2] % 22 and LinkConnected(Rd_des[-1],
                                                                                               next_node_des):  # 向后链路正常
            Rs_des.pop()
            Rd_des.append(next_node_des)
            nodes = VerticalLinksCheckForA2A(Rs_des, Rd_des, dir_hop_vertical)
            if len(nodes) > 0:  # 垂直链路正常
                new_route_from_s = Rs_des
                new_route_from_d = Rd_des
                break
        else:
            flag_2 = False  # 向后链路异常 向前结束

    return nodes, new_route_from_s, new_route_from_d


def VerticalLinksCheckForA2A(route_from_s, route_from_d, dir_hop_vertical):
    """
    垂直链路检查
    如果 链路正常 返回包括上下节点的路径，
    否则 链路异常，返回空list

    :param route_from_s:
    :param route_from_d:
    :param dir_hop_vertical:
    :return:
    """

    curr_node = route_from_s[-1]
    dest_node = route_from_d[-1]

    if curr_node // 22 != dest_node // 22:
        return []

    nodes = [curr_node]
    while curr_node != dest_node:
        next_node_m_src = Get_M_of_plane(curr_node % 22, dir_hop_vertical, '+')
        next_node = (curr_node // 22) * 22 + next_node_m_src
        if LinkConnected(curr_node, next_node):
            nodes.append(next_node)
            curr_node = next_node
        else:
            return []  # 发生错误 返回空
    return nodes


def OGDRouting_A2D(satellites, select_src_sat_id, select_des_sat_id):
    route_from_s = [select_src_sat_id]
    route_from_d = [select_des_sat_id]

    n_src = select_src_sat_id // 22
    n_des = select_des_sat_id // 22
    hop_h = n_src - n_des
    if -72 / 2 <= hop_h < 0 or hop_h > 72 / 2:
        dir_hop_horizontal = 1  # right 顺时针减 逆时针加
    elif hop_h < -72 / 2 or 0 < hop_h <= 72 / 2:
        dir_hop_horizontal = -1  # left
    else:
        dir_hop_horizontal = 0  # no movement

    m_src = select_src_sat_id % 22
    m_des = select_des_sat_id % 22
    hop_v = m_src - m_des
    if -22 / 2 <= hop_v < 0 or hop_v > 22 / 2:
        dir_hop_vertical = 1  # up 顺时针加 逆时针减
    elif hop_v < -22 / 2 or 0 < hop_v <= 22 / 2:
        dir_hop_vertical = -1  # down
    else:
        dir_hop_vertical = 0  # no movement

    while len(route_from_s) and len(route_from_d) and route_from_s[-1] % 22 != route_from_d[
        -1] % 22:  # 直到 n_src == n_des

        curr_node_src = route_from_s[-1]
        curr_node_des = route_from_d[-1]
        next_node_m_src = Get_M_of_plane(curr_node_src % 22, dir_hop_vertical, '+')
        next_node_m_des = Get_M_of_plane(curr_node_des % 22, dir_hop_vertical, '-')
        next_node_src = (curr_node_src // 22) * 22 + next_node_m_src
        next_node_des = (curr_node_des // 22) * 22 + next_node_m_des

        reward_s = abs(satellites[curr_node_src].sublat) + abs(satellites[next_node_src].sublat)
        reward_d = abs(satellites[curr_node_des].sublat) + abs(satellites[next_node_des].sublat)

        if reward_s < reward_d:
            if LinkConnected(curr_node_src, next_node_src):  # 向上链路正常
                route_from_s.append(next_node_src)
            else:  # 向上链路异常
                next_node_n_src = Get_N_of_plane(curr_node_src // 22, dir_hop_horizontal, "+")  # 向左节点
                next_node_src = next_node_n_src * 22 + curr_node_src % 22
                if LinkConnected(curr_node_src, next_node_src):  # 向左链路正常
                    route_from_s.append(next_node_src)
                else:  # 向左链路异常 回退
                    route_from_s = VerticalBacktrackingForForA2D(route_from_s, '+', dir_hop_horizontal)  # 无路可退返回空list
        else:
            if LinkConnected(curr_node_des, next_node_des):  # 正常
                route_from_d.append(next_node_des)
            else:  # 向下链路异常
                next_node_n_des = Get_N_of_plane(curr_node_des // 22, dir_hop_horizontal, '-')
                next_node_des = next_node_n_des * 22 + curr_node_des % 22
                if LinkConnected(curr_node_des, next_node_des):  # 向右链路正常
                    route_from_d.append(next_node_des)
                else:  # 向右链路异常
                    route_from_d = VerticalBacktrackingForForA2D(route_from_d, '-', dir_hop_horizontal)  # 无路可退返回空list

    if len(route_from_s) < 1 or len(route_from_d) < 1:
        return []
    elif route_from_s[-1] % 22 != route_from_d[-1] % 22:
        return []

    nodes = HorizontalLinkCheckForA2D(route_from_s, route_from_d, dir_hop_horizontal)
    if len(nodes) < 1:  # 链路故障， 链路故障 返回空为链路故障
        nodes, route_from_s, route_from_d = HorizontalSearchForA2D(route_from_s, route_from_d, dir_hop_horizontal,
                                                                   dir_hop_vertical)
        if len(route_from_s) > 0:  # 找到了路
            route_from_s += nodes[1:-1]
        else:
            return  # 无路可走
    else:  # 链路正常，路由计算完毕
        route_from_s += nodes[1:-1]

    route_from_d.reverse()
    return route_from_s + route_from_d


def VerticalBacktrackingForForA2D(route, type_from_s_or_d, dir_hop_horizontal):  # 无路可退返回空list:
    """
    向上向左链路故障时候进入，用于链路回退/回溯
    解释说明：
    for example:
    route = [1,2,3,4,5]
    temp_node->5
    curr_node->4
    5节点向上向左链路故障，回到4节点进行链路探查
    如果 temp_node 和 curr_node 具有同一轨道编号，那么从curr_node向左探查
        如果 向左探查到链路正常，将其加入路径
        否则 链路故障 再次回退
    否则 再次回退
    :param route:
    :param type_from_s_or_d:
    :param dir_hop_vertical:
    :return:
    """
    if len(route) == 1:
        return []
    temp_node = route.pop()
    curr_node = route[-1]

    if temp_node // 22 == curr_node // 22:  # 纵向退 n = n
        # 计算向左的节点编号
        next_node_n = Get_N_of_plane(curr_node // 22, dir_hop_horizontal, type_from_s_or_d)
        next_node = next_node_n * 22 + curr_node % 22
        if LinkConnected(curr_node, next_node):  # 向左正常
            route.append(next_node)
        else:  # 向左异常
            route = VerticalBacktrackingForForA2D(route, type_from_s_or_d, dir_hop_horizontal)
    else:  # 横向退，for example 4->5 说明4向上链路异常，必须再退一步
        route = VerticalBacktrackingForForA2D(route, type_from_s_or_d, dir_hop_horizontal)

    return route


def HorizontalSearchForA2D(route_from_s, route_from_d, dir_hop_horizontal, dir_hop_vertical):
    """
    在上下范围内左右寻找，
    - - - - - - - 1 d
    - - - - - 4 3 2 -
    - - - - - 5 - - -
    - - - a v b - - -
    - - 4 3 - - - - -
    - - 2 - - - - - -
    - s 2 - - - - - -
    即，在a 与 b 之间找路 采用向上向下循环找，找到一个即可
    :param route:
    :param route_from_s:
    :param route_from_d:
    :param dir_hop_horizontal:
    :param dir_hop_vertical:
    :return: new_route_from_s = []代表无路，
    """
    Rs_up = route_from_s
    Rd_up = route_from_d
    Rs_down = route_from_s
    Rd_down = route_from_d

    new_route_from_s = []
    new_route_from_d = []

    nodes = []
    flag_up = True
    flag_down = True
    while flag_up or flag_down:
        next_m_up = Get_M_of_plane(Rs_up[-1] % 22, dir_hop_vertical, '+')
        next_node_up = (Rs_up[-1] // 22) * 22 + next_m_up

        if flag_up and len(Rd_up) > 1 and Rd_up[-1] // 22 == Rd_up[-2] // 22 and LinkConnected(Rs_up[-1], next_node_up):
            Rd_up.pop()
            Rs_up.append(next_node_up)
            nodes = HorizontalLinkCheckForA2D(Rs_up, Rd_up, dir_hop_horizontal)
            if len(nodes):  # 水平链路正常
                new_route_from_s = Rs_up
                new_route_from_d = Rs_down
                break
        else:
            flag_up = False  # 向上链路异常 向上结束

        next_m_down = Get_M_of_plane(Rd_down[-1] % 22, dir_hop_vertical, '-')
        next_node_down = (Rd_down[-1] // 22) * 22 + next_m_down
        if flag_down and len(Rs_down) > 1 and Rs_down[-1] // 22 == Rs_down[-2] // 22 and LinkConnected(Rd_down,
                                                                                                       next_node_down):
            Rs_down.pop()
            Rd_down.append(next_node_down)
            nodes = HorizontalLinkCheckForA2D(Rs_down, Rd_down, dir_hop_horizontal)
            if len(nodes):  # 水平链路正常
                new_route_from_s = Rs_down
                new_route_from_d = Rd_down
                break
        else:
            flag_down = False  # 向下链路异常 向上结束

    return nodes, new_route_from_s, new_route_from_d


def HorizontalLinkCheckForA2D(route_from_s, route_from_d, dir_hop_horizontal):
    """
    水平链路检查
    如果 链路正常 返回包括左右节点的路径，
    否则 链路异常，返回空list
    :param route_from_s:
    :param route_from_d:
    :param dir_hop_horizontal:
    :return:
    """
    curr_node = route_from_s[-1]
    dest_node = route_from_d[-1]
    if curr_node % 22 != dest_node % 22:
        return []

    nodes = [curr_node]
    while curr_node != dest_node:
        next_node_n_src = Get_N_of_plane(curr_node // 22, dir_hop_horizontal, '+')
        next_node = next_node_n_src * 22 + curr_node % 22
        if LinkConnected(curr_node, next_node):
            nodes.append(next_node)
            curr_node = next_node
        else:
            return []
    return nodes


def LinkConnected(start_node, end_node) -> int:  # 链路连通

    return True


def Get_N_of_plane(n_orbit, dir, ch):
    if ch == '+':
        if n_orbit + dir < 0:
            n_orbit = 71
        else:
            n_orbit += dir
            n_orbit %= 72
        return n_orbit
    if ch == '-':
        if n_orbit - dir < 0:
            n_orbit = 71
        else:
            n_orbit -= dir
            n_orbit %= 72
        return n_orbit


def Get_M_of_plane(m_orbit, dir, ch):
    if ch == '+':
        if m_orbit + dir < 0:
            m_orbit = 21
        else:
            m_orbit += dir
            m_orbit %= 22
        return m_orbit

    if ch == '-':
        if m_orbit - dir < 0:
            m_orbit = 21
        else:
            m_orbit -= dir
            m_orbit %= 22
        return m_orbit


def orbit_gird_routing_sametype(satellites, select_src_sat_id, select_des_sat_id):
    route_from_s = [select_src_sat_id]
    route_from_d = [select_des_sat_id]

    n_src = select_src_sat_id // 22
    n_des = select_des_sat_id // 22

    m_src = select_src_sat_id % 22
    m_des = select_des_sat_id % 22

    hop_h = n_src - n_des

    if -72 / 2 <= hop_h < 0 or hop_h > 72 / 2:
        dir_hop_horizontal = 1  # right 顺时针减 逆时针加
    elif hop_h < -72 / 2 or 0 < hop_h <= 72 / 2:
        dir_hop_horizontal = -1  # left
    else:
        dir_hop_horizontal = 0  # no movement

    hop_v = m_src - m_des
    if -22 / 2 <= hop_v < 0 or hop_v > 22 / 2:
        dir_hop_vertical = 1  # up 顺时针加 逆时针减
    elif hop_v < -22 / 2 or 0 < hop_v <= 22 / 2:
        dir_hop_vertical = -1  # down
    else:
        dir_hop_vertical = 0  # no movement

    reward_s = 0
    reward_d = 0
    for i in range(abs(hop_h)):

        next_node_n_src = n_src + dir_hop_horizontal
        next_node_n_des = n_des - dir_hop_horizontal

        if n_src + dir_hop_horizontal < 0:
            next_node_n_src = 71
        if n_des - dir_hop_horizontal < 0:
            next_node_n_des = 71

        next_node_n_src %= 72
        next_node_n_des %= 72
        next_node_src = next_node_n_src * 22 + m_src
        next_node_des = next_node_n_des * 22 + m_des

        reward_s = abs(satellites[n_src * 22 + m_src].sublat) + abs(satellites[next_node_src].sublat)
        reward_d = abs(satellites[n_des * 22 + m_des].sublat) + abs(satellites[next_node_des].sublat)

        if reward_s >= reward_d:
            route_from_s.append(next_node_src)
            n_src = next_node_n_src
        else:
            route_from_d.append(next_node_des)
            n_des = next_node_n_des

    assert n_src == n_des

    for i in range(abs(hop_v) - 1):
        next_node_m_src = Get_M_of_plane(m_src, dir_hop_vertical, '+')
        next_node = n_src * 22 + next_node_m_src
        route_from_s.append(next_node)
        m_src = next_node_m_src

    route_from_d.reverse()
    return route_from_s + route_from_d


def orbit_gird_routing_difftype(satellites, select_src_sat_id, select_des_sat_id):
    route_from_s = [select_src_sat_id]
    route_from_d = [select_des_sat_id]

    n_src = select_src_sat_id // 22
    n_des = select_des_sat_id // 22

    m_src = select_src_sat_id % 22
    m_des = select_des_sat_id % 22

    hop_h = n_src - n_des

    if -72 / 2 <= hop_h < 0 or hop_h > 72 / 2:
        dir_hop_horizontal = 1  # right 顺时针减 逆时针加
    elif hop_h < -72 / 2 or 0 < hop_h <= 72 / 2:
        dir_hop_horizontal = -1  # left
    else:
        dir_hop_horizontal = 0  # no movement

    hop_v = m_src - m_des
    if -22 / 2 <= hop_v < 0 or hop_v > 22 / 2:
        dir_hop_vertical = 1  # up 顺时针加 逆时针减
    elif hop_v < -22 / 2 or 0 < hop_v <= 22 / 2:
        dir_hop_vertical = -1  # down
    else:
        dir_hop_vertical = 0  # no movement

    for i in range(hop_v):

        next_node_m_src = m_src + dir_hop_vertical
        next_node_m_des = m_des - dir_hop_vertical

        if m_src + dir_hop_vertical < 0:
            next_node_m_src = 21
        if m_des - dir_hop_vertical < 0:
            next_node_m_des = 21

        next_node_m_src %= 22
        next_node_m_des %= 22

        next_node_src = n_src * 22 + next_node_m_src
        next_node_des = n_des * 22 + next_node_m_des

        reward_s = abs(satellites[n_src * 22 + m_src].sublat) + abs(satellites[next_node_src].sublat)
        reward_d = abs(satellites[n_des * 22 + m_des].sublat) + abs(satellites[next_node_des].sublat)

        if reward_s <= reward_d:
            route_from_s.append(next_node_src)
            m_src = next_node_m_src
        else:
            route_from_d.append(next_node_des)
            m_des = next_node_m_des

    assert m_src == m_des
    route_from_d.reverse()

    for i in range(hop_h - 1):
        next_node_n_src = Get_N_of_plane(n_src, dir_hop_horizontal, '+')
        next_node = next_node_n_src * 22 + m_src
        route_from_s.append(next_node_n_src * 22 + m_src)
        n_src = next_node_n_src

    # for i in range(hop_h - 1):
    #     if n_src + dir_hop_horizontal < 0:
    #         next_node_n_src = 71
    #     else:
    #         next_node_n_src = (n_src + dir_hop_horizontal) % 72
    #     route_from_s.append(next_node_n_src * 22 + m_src)
    #     n_src = next_node_n_src

    return route_from_s + route_from_d


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
    with open(os.path.join(gen_data, sat_data[0]) + '/tles.txt', 'r') as f:
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
        "satellites": satellites
    }


def read_isls(gen_filename, num_satellites):
    """
    Read ISLs file into a list of undirected edges

    :param filename_isls:  Filename of ISLs (typically /path/to/isls.txt)
    :param num_satellites: Number of satellites (to verify indices)

    :return: List of all undirected ISL edges
    """
    isls_list = []

    info = os.listdir(gen_filename)

    with open(os.path.join(gen_filename, info[0]) + '/isls.txt', 'r') as f:
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


def ascending_or_descending_of_satellite(satellite, date_str) -> int:
    satellite.compute(date_str)
    meanAnomaly = satellite.M  # Mean anomaly from perigee at epoch
    if 90 < meanAnomaly <= 270:
        # descending satellite
        return 0
    else:
        # ascending satellite
        return 1


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
    observer.lat = str(ground_station["latitude_degrees_str"])  # Very important: string argument is in degrees.
    observer.lon = str(ground_station["longitude_degrees_str"])  # DO NOT pass a float as it is interpreted as radians
    observer.elevation = ground_station["elevation_m_float"]

    # Compute distance from satellite to observer
    satellite.compute(observer)

    # Return distance
    return satellite.range


def read_ground_stations_extended(gen_data):
    """
    Reads ground stations from the input file.

    :param filename_ground_stations_basic: Filename of ground stations basic (typically /path/to/ground_stations.txt)

    :return: List of ground stations
    """

    constellation_info_list = os.listdir(gen_data)

    ground_stations_extended = []
    gid = 0

    for i in constellation_info_list:
        with open(os.path.join(gen_data, i) + '/ground_stations.txt', 'r') as f:
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


def calculate_fstate_shortest_path_without_gs_relaying(
        satellites,
        ground_stations,
        dist_sat_net_without_gs,
        ground_station_satellites_in_range,
        enable_verbose_logs=False
):
    # Forwarding state
    fstate = {}
    # Satellites to ground stations
    # From the satellites attached to the destination ground station,
    # select the one which promises the shortest path to the destination ground station (getting there + last hop)
    dist_satellite_to_ground_station = {}
    for curr in range(len(satellites)):
        for dst_gid in range(len(ground_stations)):
            dst_gs_node_id = len(satellites) + dst_gid

            # Among the satellites in range of the destination ground station,
            # find the one which promises the shortest distance
            possible_dst_sats = ground_station_satellites_in_range[dst_gid]
            possibilities = []
            for b in possible_dst_sats:
                if not math.isinf(dist_sat_net_without_gs[(curr, b[1])]):  # Must be reachable
                    possibilities.append(
                        (
                            dist_sat_net_without_gs[(curr, b[1])] + b[0],
                            b[1]
                        )
                    )
            possibilities = list(sorted(possibilities))
            # print('possibilities', possibilities)

            # By default, if there is no satellite in range for the

            distance_to_ground_station_m = float("inf")
            if len(possibilities) > 0:
                dst_sat = possibilities[0][1]
                distance_to_ground_station_m = possibilities[0][0]
            # In any case, save the distance of the satellite to the ground station to re-use
            # when we calculate ground station to ground station forwarding
            dist_satellite_to_ground_station[(curr, dst_gs_node_id)] = distance_to_ground_station_m

        # Ground stations to ground stations
        # Choose the source satellite which promises the shortest path
    for src_gid in range(len(ground_stations)):
        for dst_gid in range(len(ground_stations)):
            if src_gid != dst_gid:
                src_gs_node_id = len(satellites) + src_gid
                dst_gs_node_id = len(satellites) + dst_gid

                # Among the satellites in range of the source ground station,
                # find the one which promises the shortest distance
                possible_src_sats = ground_station_satellites_in_range[src_gid]

                if enable_verbose_logs:
                    print('-----------------------')
                    print('----possible_src_sats----')
                    print(possible_src_sats)

                possibilities = []
                for a in possible_src_sats:
                    best_distance_offered_m = dist_satellite_to_ground_station[(a[1], dst_gs_node_id)]
                    if not math.isinf(best_distance_offered_m):
                        possibilities.append(
                            (
                                a[0] + best_distance_offered_m,
                                # distance between two stations #(51, 100): 725104.4375, (51, 101): 15010708.41926724, the two value is addde
                                a[1]
                            )
                        )
                possibilities = sorted(possibilities)

                if enable_verbose_logs:
                    print('-------------------------------')
                    print("dist_satellite_to_ground_station")
                    print(possibilities)

                # By default, if there is no satellite in range for one of the
                # ground stations, it will be dropped (indicated by -1)
                # next_hop_decision = (-1, -1, -1)
                if len(possibilities) > 0:
                    src_sat_id = possibilities[0][1]
                else:
                    print('current snapshots no gsl')
                    src_sat_id = -1

                fstate[(src_gs_node_id, dst_gs_node_id)] = src_sat_id

    return fstate


def cal_access_sats_each_time_step(ground_station_satellites_in_range) -> list:
    postions_sat_parsed = []
    for ground_station in ground_station_satellites_in_range:
        # parse num_plan and num_sat_in_plane
        # mini_dist = ground_station[0][0]
        temp = []
        for dist, id in ground_station:
            temp.append((id // NUM_SATS_PER_ORB, id % NUM_SATS_PER_ORB, dist))
        postions_sat_parsed.append(temp)

    # # def select_minimal_hop(postions_sat_parsed):
    sat_two_final_selected = []

    # minimum distance two satellites
    # sat0_orbit = 0
    # sat0_num_of_plane = 0
    # sat1_orbit = 0
    # sat1_num_of_plane = 0

    minimal_mesh_grid = 65535
    for a, b, dist0 in postions_sat_parsed[0]:
        for c, d, dist1 in postions_sat_parsed[1]:
            rows_nums = min(get_positive_int((a - c)), NUM_ORBS - get_positive_int((a - c)))
            column = min(get_positive_int((b - d)), NUM_SATS_PER_ORB - get_positive_int((b - d)))
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
                    # if dist0 + dist1 < sat0_to_gs_dist + sat1_to_gs_dist:
                    sat0_orbit = a
                    sat0_num_of_plane = b
                    final_row = rows_nums
                    sat0_to_gs_dist = dist0

                    sat1_orbit = c
                    sat1_num_of_plane = d
                    sat1_to_gs_dist = dist1

    sat_two_final_selected.append((sat0_orbit * NUM_SATS_PER_ORB + sat0_num_of_plane))
    sat_two_final_selected.append((sat1_orbit * NUM_SATS_PER_ORB + sat1_num_of_plane))
    return sat_two_final_selected


def get_positive_int(ne_or_po_num: int):
    if ne_or_po_num < 0:
        ne_or_po_num = - ne_or_po_num
    return ne_or_po_num


def main():
    print('init customized topology')
    sat_info = read_tles(gen_data)
    # print('\n>gen_data', sat_info)
    # Dictionary:{
    # "n_orbits": n_orbits,
    # "n_sats_per_orbit": n_sats_per_orbit,
    # "num_of_all_satellite": n_orbits * n_sats_per_orbit,
    # "epoch": epoch,
    # "satellites":satellites
    # }
    ground_stations = read_ground_stations_extended(gen_data)  # starlink/starlink_info/ground_stations.txt
    satellites = sat_info['satellites']
    epoch = sat_info['epoch']
    init_time = epoch + 0 * u.day
    # graph Information
    for i in range(len(satellites)):
        graphs_sat_net_graph_only_satellites_with_isls.add_node(i)
    for i in range(len(satellites) + len(ground_stations)):
        graphs_sat_net_graph_all_with_only_gsls.add_node(i)

    isl_list = read_isls(gen_data, sat_info['num_of_all_satellite'])

    for (a, b) in isl_list:
        sat_distance_m = distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(init_time))
        # if sat_distance_m > MAX_ISL_LENGTH_M:
        #     raise ValueError(
        #         "The distance between two satellites (%d and %d) "
        #         "with an ISL exceeded the maximum ISL length (%.2fm > %.2fm at t=%dns)"
        #         % (a, b, sat_distance_m, MAX_ISL_LENGTH_M, time)
        #     )

        delay = round(sat_distance_m / light_speed, 4) * 1000
        graphs_sat_net_graph_only_satellites_with_isls.add_edge(a, b, weight=delay)
        graphs_sat_net_graph_all_with_only_gsls.add_edge(a, b, weight=delay)

    if enable_verbose_logs:
        print("  > Total ISLs............. " + str(len(isl_list)))

    if enable_verbose_logs:
        print("\n  > Epoch.................. " + str(epoch))

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
            # sat_type = ascending_or_descending_of_satellite(satellites[sid], str(each_time))
            # print('type', sat_type)
            if distance_m <= MAX_GSL_LENGTH_M:
                satellites_in_range.append(
                    (distance_m, sid, ascending_or_descending_of_satellite(satellites[sid], str(epoch))))
                # graph info
                delay = round(distance_m / light_speed, 4) * 1000
                # print("\ngs2sat---->  ", delay)
                graphs_sat_net_graph_all_with_only_gsls.add_edge(sid, len(satellites) + ground_station["gid"],
                                                                 weight=delay)
        satellites_in_range = sorted(satellites_in_range)
        ground_station_satellites_in_range.append(satellites_in_range)

    if enable_verbose_logs:
        print(" \n> ground_station_satellites_in_range ")
        print(ground_station_satellites_in_range)

    routing = []
    if enable_verbose_logs:
        OGDRouting_time = []
        for i in range(0, 1):
            start_time = time.perf_counter()
            routing = OGDRouting(satellites, ground_station_satellites_in_range)
            # print('\nOGDRoute', routing)
            end_time = time.perf_counter()
            time_sum = end_time - start_time
            OGDRouting_time.append(time_sum)
        print('\nOGDRouting_time', OGDRouting_time)

    verifiy_routing(isl_list, routing[0][1], routing[0][-2], sat_info, routing[0][1:-1])

    # Calculate shortest path distances
    minWPath_vs_vt = []
    if enable_verbose_logs:
        # print("  > Calculating Floyd-Warshall for graph without ground-station relays")
        # minWPath_vs_vt = nx.dijkstra_path(self.graphs_sat_net_graph_all_with_only_gsls, source=1584, target=1585)
        # # minWPath_vs_vt_len = nx.dijkstra_path_length(mesh_net, source=source_id, target=dest_id)
        # print("\nminWPath_vs_vt", minWPath_vs_vt)

        cal_time = []
        for i in range(0, 1):
            start_time = time.perf_counter()
            # # (Note: Numpy has a deprecation warning here because of how networkx uses matrices)
            # dist_sat_net_without_gs = nx.floyd_warshall_numpy(self.graphs_sat_net_graph_all_with_only_gsls)
            minWPath_vs_vt = nx.dijkstra_path(graphs_sat_net_graph_all_with_only_gsls, source=1584,
                                              target=1585)
            # print('\nminWPath_vs_vt', minWPath_vs_vt)
            end_time = time.perf_counter()
            time_sum = end_time - start_time
            cal_time.append(time_sum)
        print("cal_time", cal_time)
    verifiy_routing(isl_list, minWPath_vs_vt[1], minWPath_vs_vt[-2], sat_info, minWPath_vs_vt[1:-1])

    # ground_fstate_dict = calculate_fstate_shortest_path_without_gs_relaying(
    # satellites,
    # ground_stations,
    # dist_sat_net_without_gs,
    # ground_station_satellites_in_range,
    # enable_verbose_logs = False)
    # print(ground_fstate_dict)


def mainTest():
    print('init customized topology')
    sat_info = read_tles(gen_data)
    # print('\n>gen_data', sat_info)
    # Dictionary:{
    # "n_orbits": n_orbits,
    # "n_sats_per_orbit": n_sats_per_orbit,
    # "num_of_all_satellite": n_orbits * n_sats_per_orbit,
    # "epoch": epoch,
    # "satellites":satellites
    # }
    ground_stations = read_ground_stations_extended(gen_data)  # starlink/starlink_info/ground_stations.txt
    satellites = sat_info['satellites']
    epoch = sat_info['epoch']
    init_time = epoch + 0 * u.day

    for i in range(len(satellites)):
        graphs_sat_net_graph_only_satellites_with_isls.add_node(i)
    for i in range(len(satellites) + len(ground_stations)):
        graphs_sat_net_graph_all_with_only_gsls.add_node(i)

    isl_list = read_isls(gen_data, sat_info['num_of_all_satellite'])

    if not os.path.exists("./output/"):
        os.makedirs("./output/")
    file_OGDRoute = open("./output/OGDRoute.txt", 'w+')
    file_OGDRouteDelay = open("./output/OGDRouteDelay.txt", 'w+')

    # total_iterations = (simulation_end_time_s / time_step_ms)
    for time_since_epoch_ms in range(0, simulation_end_time_s * 1000, time_step_ms):

        graphs_sat_net_all_node = nx.create_empty_copy(graphs_sat_net_graph_all_with_only_gsls)
        time_since_epoch = epoch + time_since_epoch_ms * u.ms
        for (a, b) in isl_list:
            sat_distance_m = distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(time_since_epoch))
            delay = round(sat_distance_m / light_speed, 4) * 1000
            graphs_sat_net_all_node.add_edge(a, b, weight=delay)

        if enable_verbose_logs:
            print("\n  > Epoch.................. " + str(epoch))
            print("\n  > Epoch.................. " + str(time_since_epoch))

        ground_station_satellites_in_range = []
        for ground_station in ground_stations:
            # Find satellites in range
            satellites_in_range = []
            for sid in range(len(satellites)):
                distance_m = distance_m_ground_station_to_satellite(
                    ground_station,
                    satellites[sid],
                    str(epoch),
                    str(time_since_epoch)
                )
                if distance_m <= MAX_GSL_LENGTH_M:
                    satellites_in_range.append(
                        (distance_m, sid, 0))
                    # graph info
                    delay = round(distance_m / light_speed, 4) * 1000
                    # print('sid', sid, 'gid', len(satellites) + ground_station["gid"])
                    graphs_sat_net_all_node.add_edge(sid, len(satellites) + ground_station["gid"], weight=delay)
            satellites_in_range = sorted(satellites_in_range)
            ground_station_satellites_in_range.append(satellites_in_range)

        if enable_verbose_logs:
            print(" \n> ground_station_satellites_in_range ")
            for gs in ground_station_satellites_in_range:
                print(gs)


        routing = []

        if enable_verbose_logs:
            print('\nOGDRoute', routing)
            temp = OGDRoutingTest(satellites, ground_station_satellites_in_range, graphs_sat_net_all_node)
            routing = [1584] + temp + [1585]

            delay = sum(graphs_sat_net_all_node[routing[i]][routing[i + 1]]['weight'] for i in range(len(routing) - 1))

            print(routing)

            for node in routing:
                file_OGDRoute.write(str(node) + ',')
            file_OGDRoute.write('\n')

            file_OGDRouteDelay.write(str(time_since_epoch_ms))
            file_OGDRouteDelay.write(",")
            file_OGDRouteDelay.write(str(delay))
            file_OGDRouteDelay.write('\n')

    # verifiy_routing(isl_list, routing[0][1], routing[0][-2], sat_info, routing[0][1:-1])
    file_OGDRoute.close()
    file_OGDRouteDelay.close()



if __name__ == "__main__":
    mainTest()
