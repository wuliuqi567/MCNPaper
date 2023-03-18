import networkx as nx
import gen_mesh_isl
import matplotlib.pyplot as plt
import random
import gen_mesh_isl


# gen_data = './gen_data2/'
# filename_ground_stations_basic_in = './gen_data2/starlink_info/ground_stations_basic.txt'
# filename_ground_stations_out = './gen_data2/starlink_info/ground_stations.txt'
# output_file = './gen_data2/starlink_info/mesh_net_path.txt'

gen_data = './test_starlink/'
filename_ground_stations_basic_in = './starlink/starlink_info/ground_stations_basic.txt'
filename_ground_stations_out = './starlink/starlink_info/ground_stations.txt'
output_file = './starlink/starlink_info/mesh_net_path.txt'

# starlink_isls_plus_grid_twostation_algorithm_free_one_only_over_isls

"""
The bandwidth of the different links is the randomly chosen from 500 KBps to 2 MBps. 
The packet loss rates of the different links are randomly chosen from 0.2 to 1.5 percent
"""


def mutipath_programming():
    """
           mesh_info = {
                 "sats_mesh_net_graph": mesh_net,
                 "delay_matrix": delay_matrix,
                 "source2dest_sats": sat_two_final_selected,
                 "mesh_net_pos": pos
                 }
    """
    mesh_info = gen_mesh_isl.get_sats_mesh_net(gen_data, filename_ground_stations_basic_in, filename_ground_stations_out)
    mesh_net = mesh_info['sats_mesh_net_graph']
    source2dest_sats = mesh_info["source2dest_sats"]
    pos = mesh_info["mesh_net_pos"]

    source_id = source2dest_sats[2]
    dest_id = source2dest_sats[3]

    sats_id = mesh_net.nodes

    write_sats_id_to_file(list(sats_id), output_file)
    # print(sats_id)
    # mesh_net.to_directed()

    minWPath_vs_vt = nx.dijkstra_path(mesh_net, source=source_id, target=dest_id)
    minWPath_vs_vt_len = nx.dijkstra_path_length(mesh_net, source=source_id, target=dest_id)
    print(minWPath_vs_vt_len)
    print(minWPath_vs_vt)
    # (gAnt.edges[edge]['weight'] for edge in nx.utils.pairwise(minWPath4)

    nx.draw(mesh_net, pos, with_labels=True, font_size=10)
    # nx.draw(mesh_net, pos, with_labels=True)
    nx.draw_networkx_nodes(mesh_net, pos, node_color='pink', nodelist=minWPath_vs_vt)

    nx.draw_networkx_nodes(mesh_net, pos, node_color='red', nodelist=[source2dest_sats[2]])
    nx.draw_networkx_nodes(mesh_net, pos, node_color='yellow', nodelist=[source2dest_sats[3]])

    # node_labels = nx.get_node_attributes(mesh_net, name='pos')
    # nx.draw_networkx_labels(mesh_net, pos, font_size=10, labels=node_labels)
    # edge_labels = nx.get_edge_attributes(mesh_net, 'bandwidth')
    # nx.draw_networkx_edge_labels(mesh_net, pos, font_size=8, edge_labels=edge_labels)

    plt.show()

#     route[[1584, 1434, 1435, 1436, 1432, 1585]]
#     > Calculating
#     Floyd - Warshall
#     for graph without ground-station relays
#
#
#     [1584, 1477, 1455, 1433, 1432, 1585]

    return 1

def get_bd_psr_weight(delay_matrix):

    # The bandwidth of the different links is the randomly chosen from 500 KBps to 2 MBps.
    # The packet loss rates of the different links are randomly chosen from 0.2 to 1.5 percent
    all_links_len = len(delay_matrix)

    bandwidth = []
    packet_lr = []
    for index in range(all_links_len):
    # bandwidth random
        band = random.randint(500, 2000)
        bandwidth.append(band)

        # packets loss rate
        plr = random.randint(2, 15) / 1000
        packet_lr.append(plr)

    return bandwidth, packet_lr




def get_qoslinks_weight(delay_matrix):

    # The bandwidth of the different links is the randomly chosen from 500 KBps to 2 MBps.
    # The packet loss rates of the different links are randomly chosen from 0.2 to 1.5 percent

    qos = [[0.33, 0.33, 0.33],
           [0.1, 0.8, 0.1],
           [0.8, 0.1, 0.1]]

    all_links_len = len(delay_matrix)

    bandwidth = []
    packet_lr = []
    for index in range(all_links_len):
    # bandwidth random
        band = random.randint(500, 2000)
        bandwidth.append(band)

        # packets loss rate
        plr = random.randint(2, 15) / 1000
        packet_lr.append(plr)

    max_delay = max(delay_matrix)
    min_delay = min(delay_matrix)

    max_bandwidth = max(bandwidth)
    min_bandwidth = min(bandwidth)

    max_packet_lr = max(packet_lr)
    min_packet_lr = min(packet_lr)

    for index in range(all_links_len):
        delay_matrix[index] = (delay_matrix[index] - min_delay) / (max_delay - min_delay)
        packet_lr[index] = (packet_lr[index] - min_packet_lr) / (max_packet_lr - min_packet_lr)
        bandwidth[index] = (bandwidth[index] - min_bandwidth) / (max_bandwidth - min_bandwidth)

    Qb = round(min_bandwidth / sum(bandwidth), 4)
    Qd = round(min_delay / sum(delay_matrix), 4)
    Ql = round(min_packet_lr / sum(packet_lr), 6)

    qos_links_weight =[]

    qos0link_weith = []
    qos0 = qos[0]
    for index in range(all_links_len):
        weight = round(qos0[0]*Qb*bandwidth[index] - qos0[1]*Qd*delay_matrix[index] - qos0[2]*Ql*packet_lr[index], 2)
        if weight<=0:
            weight = 0.01
        qos0link_weith.append(weight)
    qos_links_weight.append(maxminnormalization(qos0link_weith))


    qos1link_weith = []
    qos1 = qos[1]
    for index in range(all_links_len):
        weight = round(qos1[0]*Qb*bandwidth[index] - qos1[1]*Qd*delay_matrix[index] - qos1[2]*Ql*packet_lr[index], 2)
        if weight <= 0:
            weight = 0.01
        qos1link_weith.append(weight)

    qos_links_weight.append(maxminnormalization(qos1link_weith))


    qos2link_weith = []
    qos2 = qos[2]
    for index in range(all_links_len):
        weight = round(qos2[0]*Qb*bandwidth[index] - qos2[1]*Qd*delay_matrix[index] - qos2[2]*Ql*packet_lr[index], 2)
        if weight <= 0:
            weight = 0.01
        qos2link_weith.append(weight)
    qos_links_weight.append(maxminnormalization(qos2link_weith))

    return qos_links_weight



def write_sats_id_to_file(sats_list, output_file):
    with open(output_file, 'w') as f:
        for node in sats_list:
            f.write(str(node))
            f.write(',')




if __name__=='__main__':
    mutipath_programming()

    # gs_stations_file = './gen_data2/starlink_info/gs_with_time.txt'
    # gen_mesh_isl.get_sats_over_gs_for_time(gs_stations_file)

#final-path-new-2-los
# 9,New-York-Newark,40.717042,-74.003663,0
# 20,Los-Angeles-Long-Beach-Santa-Ana,34.031656,-118.241716,0
# two_citys = [9, 20]
# best optimal path
# 1477,1455,1433,1411,1389,1367
# 1477,1476,1454,1432,1410,1388,1366,1367
# 1477,1478,1456,1434,1412,1390,1368,1367
#mesh_net_path.txt
# 1366,1367,1368,1388,1389,1390,1410,1411,1412,1432,1433,1434,1454,1455,1456,1476,1477,1478,

#d  algorithm
# path
# 1085,1084,1083,1105,1127,1149,1171,1193,1215,1237,1259,1281,1303,1325,1347,1369,1368,1367
# mesh_net_path.txt
# 1058,1059,1060,1061,1062,1063,1064,1080,1081,1082,1083,1084,1085,1086,1102,1103,1104,1105,1106,1107,1108,1124,1125,1126,1127,1128,1129,1130,1146,1147,1148,1149,1150,1151,1152,1168,1169,1170,1171,1172,1173,1174,1190,1191,1192,1193,1194,1195,1196,1212,1213,1214,1215,1216,1217,1218,1234,1235,1236,1237,1238,1239,1240,1256,1257,1258,1259,1260,1261,1262,1278,1279,1280,1281,1282,1283,1284,1300,1301,1302,1303,1304,1305,1306,1322,1323,1324,1325,1326,1327,1328,1344,1345,1346,1347,1348,1349,1350,1366,1367,1368,1369,1370,1371,1372,1388,1389,1390,1391,1392,1393,1394,
# delay = 67.6ms

# 123, San - Francisco - Oakland, 37.759881, -122.437392, 0
# 32, Bogotá, 4.60971, -74.08175, 0
# two_citys = [123, 32]

# 844,866,888,910,911,912
# 844,845,867,889,890,912
# 844,843,865,887,909,910,911,912
# 844,822,823,824,846,868,890,912
# d algorithm
#1323,1324,1325,1303,1281,1259,1237,1215,1193,1171,1149,1127,1105,1083,1061,1039,1017,995,973,951,952,953,954,955,956












    # path_plr = [mesh_net.edges[edge]['plr'] for edge in nx.utils.pairwise(minWPath_vs_vt)]
    # path_bd = [mesh_net.edges[edge]['bandwidth'] for edge in nx.utils.pairwise(minWPath_vs_vt)]
    # delay_weight = [mesh_net.edges[edge]['weight'] for edge in nx.utils.pairwise(minWPath_vs_vt)]
    #
    # path_plr2 = gen_mesh_isl.maxminnormalization(path_plr)
    # path_bd2 = gen_mesh_isl.maxminnormalization(path_bd)
    # delay_weight2 = gen_mesh_isl.maxminnormalization(delay_weight)
    #
    # qos = [[0.33, 0.33, 0.33],
    #        [0.1, 0.8, 0.1],
    #        [0.8, 0.1, 0.1]]
    # quilty = []
    # for i in range(len(qos)):
    #     temp = qos[i][0] * sum(path_bd2) - qos[i][1] * sum(delay_weight2) - qos[i][2] * sum(path_plr2)
    #     quilty.append(temp)
    #
    # min2path = [844, 845, 867, 868, 869, 891]
    #
    # path_plr1 = [mesh_net.edges[edge]['plr'] for edge in nx.utils.pairwise(min2path)]
    # path_bd1 = [mesh_net.edges[edge]['bandwidth'] for edge in nx.utils.pairwise(min2path)]
    # delay_weight1 = [mesh_net.edges[edge]['weight'] for edge in nx.utils.pairwise(min2path)]
    #
    # path_plr11 = gen_mesh_isl.maxminnormalization(path_plr1)
    # path_bd11 = gen_mesh_isl.maxminnormalization(path_bd1)
    # delay_weight11 = gen_mesh_isl.maxminnormalization(delay_weight1)
    #
    # for i in range(len(qos)):
    #     temp = qos[i][0] * sum(path_bd11) - qos[i][1] * sum(delay_weight11) - qos[i][2] * sum(path_plr11)
    #     quilty.append(temp)
    #
    # min3path = [844, 822, 823, 824, 825, 847, 869, 891]
    #
    # path_plr3 = [mesh_net.edges[edge]['plr'] for edge in nx.utils.pairwise(min3path)]
    # path_bd3 = [mesh_net.edges[edge]['bandwidth'] for edge in nx.utils.pairwise(min3path)]
    # delay_weight3 = [mesh_net.edges[edge]['weight'] for edge in nx.utils.pairwise(min3path)]
    #
    # path_plr33 = gen_mesh_isl.maxminnormalization(path_plr3)
    # path_bd33 = gen_mesh_isl.maxminnormalization(path_bd3)
    # delay_weight33 = gen_mesh_isl.maxminnormalization(delay_weight3)
    #
    # for i in range(len(qos)):
    #     temp = qos[i][0] * sum(path_bd33) - qos[i][1] * sum(delay_weight33) - qos[i][2] * sum(path_plr33)
    #     quilty.append(temp)
    #
    # lMinWPath1 = nx.dijkstra_path_length(mesh_net, source=source_id, target=dest_id)  # 最短加权路径长度
    #
    # all_path = nx.all_simple_paths(mesh_net, source=source_id, target=dest_id)
    # print(minWPath_vs_vt)
    # print(lMinWPath1)
    # for path in all_path:
    #     print(path)

