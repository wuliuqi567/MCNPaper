# networkX_E3.py
# Demo of shortest path with NetworkX
# Copyright 2021 YouCans, XUPT
# Crated：2021-05-20
#
# import matplotlib.pyplot as plt # 导入 Matplotlib 工具包
# import networkx as nx  # 导入 NetworkX 工具包
#
# # 问题 1：蚂蚁的最优路径分析（西安邮电大学第12届数学建模竞赛B题）
#
# gAnt = nx.Graph()  # 创建：空的 无向图
# gAnt.add_weighted_edges_from([(0,1,3),(0,2,1),(0,3,1),
#                             (1,2,1),(1,4,1),(1,9,4),
#                             (2,3,1),(2,4,2),(2,5,1),
#                             (3,5,2),(3,6,2),(3,7,1),
#                             (4,5,1),(4,9,1),
#                             (5,6,1),(5,9,3),(5,10,1),(5,12,3),
#                             (6,7,1),(6,8,2),(6,12,2),(6,13,4),(6,14,3),
#                             (7,8,1),
#                             (8,14,1),(8,15,3),
#                             (9,10,1),(9,11,1),
#                             (10,11,1),(10,12,2),
#                             (11,12,1),(11,16,1),
#                             (12,13,2),(12,16,1),
#                             (13,14,1),(13,15,2),(13,16,2),(13,17,1),
#                             (14,15,1),
#                             (15,17,4),
#                             (16,17,1)])  # 向图中添加多条赋权边: (node1,node2,weight)
#
# pos={0:(1,8),1:(4,12),2:(4,9),3:(4,6),4:(8,11),5:(9,8),  # 指定顶点位置
#      6:(11,6),7:(8,4),8:(12,2),9:(12,13),10:(15,11),11:(18,13),
#      12:(19,9),13:(22,6),14:(18,4),15:(21,2),16:(22,11),17:(28,8)}
# nx.draw(gAnt, pos, with_labels=True, alpha=0.8)
# labels = nx.get_edge_attributes(gAnt,'weight')
# nx.draw_networkx_edge_labels(gAnt,pos,edge_labels=labels, font_color='c') # 显示权值
# nx.draw_networkx_nodes(gAnt,pos,nodelist=[0,17],node_color='yellow')  # 设置顶点颜色
# nx.draw_networkx_nodes(gAnt,pos,nodelist=[7,12],node_color='lime')  # 设置顶点颜色
# nx.draw_networkx_edges(gAnt,pos,edgelist=[(2,4),(13,14)],edge_color='lime',width=2.5)  # 设置边的颜色
# nx.draw_networkx_edges(gAnt,pos,edgelist=[(11,12)],edge_color='r',width=2.5)  # 设置边的颜色
# plt.show()
#
#
# al_path = [path for path in nx.all_simple_paths(gAnt, 0, 17)]
# print(al_path)
#
# # 4. 限制条件：多个必经点 (N7,N15)
# # 解决方案：遍历从起点到终点的简单路径，求满足必经点条件的最短路径
# minWPath4 = min([path  # 返回 key 为最小值的 path
#     for path in nx.all_simple_paths(gAnt, 0, 17)  # gAnt 中所有起点为0、终点为17的简单路径
#     if all(n in path for n in (7, 15))], # 满足路径中包括顶点 N7,N15
#     key=lambda x: sum(gAnt.edges[edge]['weight'] for edge in nx.utils.pairwise(x))) # key 为加权路径长度
# lenPath = sum(gAnt.edges[edge]['weight'] for edge in nx.utils.pairwise(minWPath4))  # 求指定路径的加权路径长度
# print("\n问题4: 多个必经点的约束")
# print("S 到 E 的最短加权路径: ", minWPath4)
# print("S 到 E 的最短加权路径长度: ", lenPath)



