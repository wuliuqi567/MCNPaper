import sys
sys.path.append("../satgenpy")
import satgen
import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import openpyxl




directory = './Globalstar'
# 打开 Excel 文件
workbook = openpyxl.load_workbook('access.xlsx')  # 替换为你的文件路径

# 选择工作表（这里假设选择第一个工作表）
sheet_iridium = workbook['iridium']
sheet_Globalstar = workbook['Globalstar']
sheet_Glonass = workbook['Glonass']


# 使用os.listdir()函数获取目录下所有文件和子目录的列表
file_list = os.listdir(os.path.join(directory, 'isls'))
file_list = sorted(file_list)
# 遍历文件列表，过滤出文件并打印

service_weight = [[0.4, 0.5, 0.1],
                  [0.65, 0.15, 0.2],
                  [0.25, 0.45, 0.3]]
# 行控制 每一行代表一个时间拓扑
row_idx = 1
for file in file_list:
    walker_star_network = nx.Graph()
    with open(os.path.join(directory, 'isls', file), 'r') as f:
        lines = f.readlines()

        for line in lines[1::]:
            line = line.strip().split()
            src_node = int(line[0])
            des_node = int(line[1])
            band = float(line[2])
            delay = float(line[3])
            plr = float(line[4])
            walker_star_network.add_edge(src_node, des_node, weight=delay, band=band, plr=plr)

    conn = list(sheet_Globalstar.iter_rows())[row_idx]
    row_idx += 1
    service_1 = {'src': conn[0].value, 'dest': conn[1].value, 'bandwidth': 100/8, 'delay': 0.5, 'loss_rate': 0.3}
    service_2 = {'src': conn[2].value, 'dest': conn[3].value, 'bandwidth': 700/8, 'delay': 0.7, 'loss_rate': 0.2}
    service_3 = {'src': conn[4].value, 'dest': conn[5].value, 'bandwidth': 400/8, 'delay': 0.15, 'loss_rate': 0.15}

    s1_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_1['src'], service_1['dest'])
    s2_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_2['src'], service_2['dest'])
    s3_shortestLenghtPath = nx.dijkstra_path(walker_star_network, service_3['src'], service_3['dest'])
    # print("shortest")
    # print(s1_shortestLenghtPath)
    # print(s2_shortestLenghtPath)
    # print(s3_shortestLenghtPath)
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

    max_u1 = -1000
    max_u2 = -1000
    max_u3 = -1000
    path_max_u1 = []
    path_max_u2 = []
    path_max_u3 = []
    # print(paths)

    def osr(walker_star_network, one_path, rank):
        # print(one_path)
        path_list = []
        for i in range(len(one_path) - 1):
            path_list.append(sorted([one_path[i], one_path[i + 1]]))

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

        # print(str(delay))
        # print(max_d, min_d)
        for i in range(len(delay)):
            if max_d != min_d:
                delay[i] = round((delay[i] - min_d) / (max_d - min_d), 6)
            else:
                delay[i] = 0.5
            if max_b != min_b:
                band[i] = round((band[i] - min_b) / (max_b - min_b), 6)
            else:
                band[i] = 0.5
            # print(max_p, min_p)
            if max_p != min_p:
                plr[i] = round((plr[i] - min_p) / (max_p - min_p), 6)
            else:
                plr[i] = 0.5

        # min_b = min(band)
        # max_d = max(delay)
        # max_p = max(plr)

        u = service_weight[rank][0] * qb * sum(band) - service_weight[rank][1] * qd * sum(delay) - service_weight[rank][
            2] * qp * sum(plr)
        return u

    for path in list(s1_paths):
        # print(path)

        # print('one', one_path)
        u1 = osr(walker_star_network, path, 0)
        # print(u)
        if u1 > max_u1:
            max_u1 = u1
            path_max_u1 = path

    for path in list(s2_paths):
        # print(path)

        # print('one', one_path)
        u2 = osr(walker_star_network, path, 1)
        # print(u)
        if u2 > max_u2:
            max_u2 = u2
            path_max_u2 = path

    for path in list(s3_paths):
        # print(path)

        # print('one', one_path)
        u3 = osr(walker_star_network, path, 2)
        # print(u)
        if u3 > max_u3:
            max_u3 = u3
            path_max_u3 = path

    with open(os.path.join(directory, 'path.txt'), 'a') as f:
        f.write(str(path_max_u1))
        f.write('\n')
        f.write(str(path_max_u2))
        f.write('\n')
        f.write(str(path_max_u3))
        f.write('\n')











#     print(nx.path_weight(walker_network, path))
#     length = nx.shortest_path_length(walker_network, service_1['src'], service_1['dest'])
#     print(length)

# for path in map(nx.utils.pairwise, paths):
#     one_path = list(path)



# nx.draw(walker_network, with_labels=True, font_weight='bold')
# edge_labels = nx.get_edge_attributes(walker_network, 'bandwidth')
# # nx.draw_networkx_edge_labels(walker_network, edge_labels=edge_labels)
# plt.show()





