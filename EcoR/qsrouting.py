import sys
sys.path.append("../satgenpy")
import satgen
import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


# # read configure file
# configure_dict = {}
# with open('config.txt', 'r') as config_file:
#     lines = config_file.readlines()
#     for line in lines:
#         line = line.strip().split('=')
#         configure_dict.setdefault(line[0], line[1])
#
# gen_data = configure_dict['satellite_network_dir_and_name']
#
# tle_file = os.path.join(gen_data, "tles.txt")
# isl_file = os.path.join(gen_data, "isls.txt")
# distance_file = os.path.join(gen_data, "distance.txt")
#
# sat_info = satgen.read_tles(tle_file)
# satellites = sat_info['satellites']
# epoch = sat_info['epoch']
#
# isl_list = satgen.read_isls(isl_file, len(satellites))
#
# walker_network = nx.Graph()
# for i in range(len(satellites)):
#     walker_network.add_node(i)
#
# for (a, b) in isl_list:
#     sat_distance_m = satgen.distance_m_between_satellites(satellites[a], satellites[b], str(epoch), str(epoch))
#     delay = round(sat_distance_m/3e8*1000, 1)
#     # print(delay)
#     walker_network.add_edge(a, b, weight=delay, band=np.random.randint(100, 110), plr=round(0.0015 + 0.0005*np.random.rand(), 8))


# indices0 = "9 13 26 6 5 7 8 14 1 21 24 15 20 19 4 25 2 27 3 23 11 22 12 28 30 29 16 17 10 18"
# indices1 = "2 9 13 19 16 12 14 24 4 29 8 25 30 20 28 17 1 22 7 21 11 23 6 3 26 5 15 18 27 10"
# indices2 = "26 3 23 29 24 11 16 1 18 9 13 6 30 21 7 4 2 27 19 15 20 28 10 14 25 12 22 8 17 5"
# indices3 = "1 19 25 7 17 16 6 11 30 4 9 21 24 28 8 15 29 13 12 14 2 27 18 5 26 20 22 3 23 10"
# indices4 = "16 12 30 24 15 4 26 8 5 10 6 11 18 2 1 20 9 27 28 19 23 14 29 25 21 3 17 13 22 7"
# indices5 = "19 24 18 22 12 26 6 25 3 10 20 14 21 16 9 5 28 23 11 30 17 8 29 4 1 13 15 7 27 2"
#
# indices0 = indices0.split(" ")
# indices1 = indices1.split(" ")
# indices2 = indices2.split(" ")
# indices3 = indices3.split(" ")
# indices4 = indices4.split(" ")
# indices5 = indices5.split(" ")

walker_star_network = nx.Graph()
for i in range(1, 67):
    walker_star_network.add_node(i)


with open("./tu2.csv", 'r') as f:
    lines = f.readlines()

    for line in lines:
        line = line.strip().split(',')
        src_node = int(line[0])
        des_node = int(line[1])
        delay = float(line[3])
        band = float(line[4])
        plr = float(line[5])

        walker_star_network.add_edge(src_node, des_node, weight=delay, band=band, plr=plr)


service_1 = {'src': 11, 'dest': 56, 'bandwidth': 100/8, 'delay': 0.5, 'loss_rate': 0.3}
service_2 = {'src': 22, 'dest': 49, 'bandwidth': 700/8, 'delay': 0.7, 'loss_rate': 0.2}
service_3 = {'src': 33, 'dest': 55, 'bandwidth': 400/8, 'delay': 0.15, 'loss_rate': 0.15}

service_weight = [[0.4, 0.5, 0.1],
                  [0.65, 0.15, 0.2],
                  [0.25, 0.45, 0.3]]


s1_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_1['src'], service_1['dest'])
s2_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_2['src'], service_2['dest'])
s3_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_2['src'], service_2['dest'])
print("shortest")
print(s1_shortestLenghtPath)
print(s2_shortestLenghtPath)
print(s3_shortestLenghtPath)
# shortestLenght = nx.shortest_path_length(walker_network, service_1['src'], service_1['dest'])
# print(shortestLenghtPath, shortestLenght)

s1_minLengthPath = len(s1_shortestLenghtPath) - 1
s1_maxLengthPath = len(s1_shortestLenghtPath) + 1

s2_minLengthPath = len(s2_shortestLenghtPath) - 1
s2_maxLengthPath = len(s2_shortestLenghtPath) + 1

s3_minLengthPath = len(s3_shortestLenghtPath) - 1
s3_maxLengthPath = len(s3_shortestLenghtPath) + 1


# 而cutoff参数设置的截断长度不包括源结点。
s1_paths = nx.all_simple_paths(walker_star_network, service_1['src'], service_1['dest'], cutoff=s1_minLengthPath)
s1_paths = list(s1_paths)
s2_paths = nx.all_simple_paths(walker_star_network, service_2['src'], service_2['dest'], cutoff=s2_minLengthPath)
s2_paths = list(s2_paths)
s3_paths = nx.all_simple_paths(walker_star_network, service_3['src'], service_3['dest'], cutoff=s3_minLengthPath)
s3_paths = list(s3_paths)


def osr(one_path, rank):

    path_list = []
    for i in range(len(one_path)-1):
        path_list.append(sorted([one_path[i], one_path[i+1]]))

    # print(one_path)
    # print(path_list)

    delay = []
    band = []
    plr = []
    all_delay = nx.get_edge_attributes(walker_star_network, "weight")
    all_band = nx.get_edge_attributes(walker_star_network, "band")
    all_plr = nx.get_edge_attributes(walker_star_network, "plr")

    for i in range(len(path_list)):
        # aa = all_delay.get(tuple(path_list[i]))
        if all_delay.get(tuple(path_list[i])):
            delay.append(all_delay[tuple(path_list[i])])
            band.append(all_band[tuple(path_list[i])])
            plr.append(all_plr[tuple(path_list[i])])
        else:
            # print('hh')
            link = path_list[i]
            link.reverse()
            delay.append(all_delay[tuple(link)])
            band.append(all_band[tuple(link)])
            plr.append(all_plr[tuple(link)])

    # normal
    min_d = min(delay)
    min_b = min(band)
    min_p = min(plr)

    max_d = max(delay)
    max_b = max(band)
    max_p = max(plr)

    qd = round(max(delay) / sum(delay), 6)
    qb = round(min(band) / sum(band), 6)
    qp = round(max(plr) / sum(plr), 6)

    for i in range(len(delay)):
        delay[i] = round((delay[i] - min_d) / (max_d - min_d), 6)
        band[i] = round((band[i] - min_b) / (max_b - min_b), 6)
        # print(max_p, min_p)
        plr[i] = round((plr[i] - min_p) / (max_p - min_p), 6)

    # min_b = min(band)
    # max_d = max(delay)
    # max_p = max(plr)

    u = service_weight[rank][0] * qb * sum(band) - service_weight[rank][1] * qd * sum(delay) - service_weight[rank][2] * qp * sum(plr)
    return u


max_u1 = -1000
max_u2 = -1000
max_u3 = -1000
path_max_u1 =[]
path_max_u2 =[]
path_max_u3 =[]
# print(paths)


for path in list(s1_paths):
    # print(path)

    # print('one', one_path)
    u = osr(path, 0)
    # print(u)
    if u > max_u1:
        max_u1 = u
        path_max_u1 = path

for path in list(s2_paths):
    # print(path)

    # print('one', one_path)
    u = osr(path, 0)
    # print(u)
    if u > max_u1:
        max_u1 = u
        path_max_u2 = path

for path in list(s3_paths):
    # print(path)

    # print('one', one_path)
    u = osr(path, 0)
    # print(u)
    if u > max_u1:
        max_u1 = u
        path_max_u3 = path

print("osr")
print(path_max_u1)
print(path_max_u2)
print(path_max_u3)





#     print(nx.path_weight(walker_network, path))
#     length = nx.shortest_path_length(walker_network, service_1['src'], service_1['dest'])
#     print(length)

# for path in map(nx.utils.pairwise, paths):
#     one_path = list(path)



# nx.draw(walker_network, with_labels=True, font_weight='bold')
# edge_labels = nx.get_edge_attributes(walker_network, 'bandwidth')
# # nx.draw_networkx_edge_labels(walker_network, edge_labels=edge_labels)
# plt.show()





