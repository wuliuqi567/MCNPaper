### process

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

### 配置信息

starlink

```
satellite_network_dir_and_name=test_starlink
ALTITUDE_M=550000
Elevation_angle=25
NUM_ORBS=72
NUM_SATS_PER_ORB=22
INCLINATION_DEGREE=53
simulation_end_time_s=600
Time_step_ms=100
isl_data_rate_megabit_per_s=10.0
gsl_data_rate_megabit_per_s=10.0
isl_max_queue_size_pkts=100
gsl_max_queue_size_pkts=100
```

```
satellite_network_dir_and_name=test_kuiper
ALTITUDE_M=630000
Elevation_angle=30
NUM_ORBS=34
NUM_SATS_PER_ORB=34
INCLINATION_DEGREE=52
simulation_end_time_s=600
Time_step_ms=100
isl_data_rate_megabit_per_s=10.0
gsl_data_rate_megabit_per_s=10.0
isl_max_queue_size_pkts=100
gsl_max_queue_size_pkts=100
```

```
satellite_network_dir_and_name=test_telesat
ALTITUDE_M=630000
Elevation_angle=10
NUM_ORBS=27
NUM_SATS_PER_ORB=13
INCLINATION_DEGREE=99
simulation_end_time_s=600
Time_step_ms=100
isl_data_rate_megabit_per_s=10.0
gsl_data_rate_megabit_per_s=10.0
isl_max_queue_size_pkts=100
gsl_max_queue_size_pkts=100

```


## Getting started

1. Python 3.7+
2. The following dependencies need to be installed:
   ```
   pip install numpy astropy ephem networkx sgp4 geopy matplotlib statsmodels
   sudo apt-get install libproj-dev proj-data proj-bin libgeos-dev
   pip install cartopy
   pip install git+https://github.com/snkas/exputilpy.git@v1.6
   ```

zenme buxing
