# -*- coding:utf8 -*-
import os
import time
import numpy as np
import scipy.sparse.linalg
from scipy import linalg
# ![](../../../../../var/folders/zs/s5ymqv5913z46y8s93shrt_00000gn/T/cc.ffitch.shottr/cc.ffitch.shottr/SCR-20230501-siul.png)

def main():
    start_time = time.time()


    ######## 矩阵求解的b start
    every_people_on_plant = "./people_on_plant.txt"
    f2 = open(every_people_on_plant)
    lines2 = f2.readlines()
    matrix_b = np.ones(len(lines2) * 2)
    dic2 = {}
    for i, line in enumerate(lines2):
        dic2[line.split(":")[0]] = i
        matrix_b[i] = line.split(":")[1].split("\n")[0].split(",")[0]
        matrix_b[i + len(lines2)] = line.split(":")[1].split("\n")[0].split(",")[1]
    ##### 矩阵求解的b end

    ######## 矩阵求解的a start
    every_path = "./path.txt"
    f1 = open(every_path)
    lines1 = f1.readlines()
    matrix_a = np.zeros([2 * len(lines2), len(lines1)])
    for i, line in enumerate(lines1):
        for j, s in enumerate(line.split(":")[1].split("\n")[0].split(",")):
            matrix_a[(j % 2) * len(lines2) + dic2[s]][i] = 1
    ##### 矩阵求解的a end
    # print(matrix_a[145][7])
    # print(matrix_a[148 + 1883][7])
    # print(matrix_b[145])
    # print(matrix_b[148 + 1883])
    ### 矩阵求解
    x = scipy.sparse.linalg.lsmr(matrix_a, matrix_b)[0]
    y = []
    for k in range(len(x)):
        y.append("path_{}:{}\n".format(k + 1, x[k]))

    with open(data_path + one_data + "/result.txt", "w") as f:
        f.writelines(y)
    # dic = {plant_entry:i for i, plant_entry in enumerate()}
    # process();


if __name__ == "__main__":
    main()
