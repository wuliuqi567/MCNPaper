


###  process
1. generate mini-net satellites---> lino_mininet_topo.py
2. add all flow ----_> add_flow.py
3. config_mininet.properties configure file

说明

卫星tle数据生成 main.py  main_helper.py  gen_satellies_info.py
配置文件 config_mininet.csv

路由计算 starlinkRouting.py

画图  mesh_mutipath_programming.py gen_mesh_isl



判断上升与下降//

选择跳数最少的接入卫星

链路只有通断两种属性，执行路径计算



reference：
Analysis of Inter-Satellite Link Paths for LEO Mega-Constellation Networks
Distributed On-Demand Routing for LEO Mega-Constellations: A Starlink Case Study

Dynamic Routing for Software-Defined LEO Satellite Networks based on ISL Attributes