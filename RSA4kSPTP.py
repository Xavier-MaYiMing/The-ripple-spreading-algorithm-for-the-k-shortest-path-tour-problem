#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/1/2 15:47
# @Author  : Xavier Ma
# @Email   : xavier_mayiming@163.com
# @File    : RSA4kSPTP.py
# @Statement : The ripple-spreading algorithm for the k shortest path tour problem
import copy
import numpy as np
import matplotlib.pyplot as plt
import random
import math


def find_neighbor(network):
    """
    Find the neighbor of each node
    :param network:
    :return: [[the neighbor nodes of node 1], ...]
    """
    nn = len(network)
    neighbor = []
    for i in range(nn):
        neighbor.append(list(network[i].keys()))
    return neighbor


def find_speed(network, neighbor):
    """
    Find the ripple-spreading speed
    :param network:
    :param neighbor:
    :return:
    """
    speed = 1e10
    for i in range(len(network)):
        for j in neighbor[i]:
            speed = min(speed, network[i][j])
    return speed


def cal_cost(network, path):
    """
    calculate the cost of the path
    :param network:
    :param path:
    :return:
    """
    cost = 0
    for i in range(len(path) - 1):
        cost += network[path[i]][path[i + 1]]
    return cost


def subRSA(network, neighbor, source, destination, init_time, init_radius, init_path, v, k):
    """
    The ripple-spreading algorithm for the subproblems of the k-SPTP
    :param network: {node1: {node2: length, node3: length, ...}, ...}
    :param neighbor: the neighbor set
    :param source: the set of source nodes
    :param destination: the set of destination nodes
    :param init_time: the initial time for each initial ripple
    :param init_radius: the initial radius for each initial ripple
    :param init_path: the initial path for each initial ripple
    :param v: the ripple-spreading speed
    :param k: the k shortest paths
    :return:
    """
    # Step 1. Initialization
    nn = len(network)  # node number
    t = min(init_time) - 1
    nr = 0  # the number of ripples - 1
    epicenter_set = []  # epicenter set
    radius_set = []  # radius set
    path_set = []  # path set
    active_set = []  # the set containing all active ripples
    start_flag = source.copy()
    dest_ripple = {}  # the ripple reaching destinations
    for node in destination:
        dest_ripple[node] = []
    omega = {}  # the set that records the ripple generated at each node
    for node in range(nn):
        omega[node] = []

    # Step 2. The main loop
    while True:

        # Step 2.1. Termination judgment
        flag = True
        for node in destination:
            if len(omega[node]) < k:  # there is a destination node that has not been visited k times by ripples
                flag = False
                break
        if flag:
            break

        # Step 2.2. Time updates
        t += 1
        incoming_ripples = {}
        new_path = []
        for ripple in active_set:

            # Step 2.3. Active ripples spread out
            radius_set[ripple] += v

            # Step 2.4. New incoming ripples
            radius = radius_set[ripple]
            epicenter = epicenter_set[ripple]
            path = path_set[ripple]
            for node in neighbor[epicenter]:
                if len(omega[node]) < k:  # the node is visited no more than k times
                    temp_length = network[epicenter][node]
                    if temp_length <= radius < temp_length + v:
                        temp_path = path.copy()
                        temp_path.append(node)
                        if temp_path not in new_path:
                            new_path.append(temp_path)
                            if node in incoming_ripples.keys():
                                incoming_ripples[node].append({
                                    'path': temp_path,
                                    'radius': radius - temp_length,
                                })
                            else:
                                incoming_ripples[node] = [{
                                    'path': temp_path,
                                    'radius': radius - temp_length,
                                }]

        # Step 2.5. Generate initial ripples
        if start_flag:
            need_to_delete = []
            for i in range(len(start_flag)):
                if t == init_time[i]:
                    need_to_delete.append(i)
                    node = start_flag[i]
                    if len(omega[node]) < k and init_path[i] not in new_path:
                        new_path.append(init_path[i])
                        if node in incoming_ripples.keys():
                            incoming_ripples[node].append({
                                'path': init_path[i],
                                'radius': init_radius[i],
                            })
                        else:
                            incoming_ripples[node] = [{
                                'path': init_path[i],
                                'radius': init_radius[i],
                            }]
            for i in range(len(need_to_delete) - 1, -1, -1):
                start_flag.pop(i)
                init_time.pop(i)
                init_radius.pop(i)
                init_path.pop(i)

        # Step 2.6. Trigger new ripples
        for node in incoming_ripples.keys():
            new_ripples = sorted(incoming_ripples[node], key=lambda x: x['radius'], reverse=True)
            if len(omega[node]) + len(new_ripples) > k:
                new_ripples = new_ripples[: k - len(omega[node])]
            for item in new_ripples:
                path_set.append(item['path'])
                epicenter_set.append(node)
                radius_set.append(item['radius'])
                active_set.append(nr)
                omega[node].append(nr)
                nr += 1
                if node in destination:
                    dest_ripple[node].append({
                        'radius': item['radius'],
                        'time': t,
                        'path': item['path'],
                    })

        # Step 2.7. Active -> Inactive
        remove_ripple = []
        for ripple in active_set:
            epicenter = epicenter_set[ripple]
            radius = radius_set[ripple]
            flag_inactive = True
            for node in neighbor[epicenter]:
                if radius < network[epicenter][node] and len(omega[node]) < k:
                    flag_inactive = False
                    break
            if flag_inactive:
                remove_ripple.append(ripple)
        for ripple in remove_ripple:
            active_set.remove(ripple)

    # Step 3. Sort the results
    dest_node = []
    dest_time = []
    dest_radius = []
    dest_path = []
    for node in destination:
        for item in dest_ripple[node]:
            dest_node.append(node)
            dest_path.append(item['path'])
            dest_time.append(item['time'])
            dest_radius.append(item['radius'])
    return dest_node, dest_time, dest_radius, dest_path


def main(network, node_subset, k):
    """
    The main function of the RSA4kSPTP
    :param network: {node1: {node2: length, node3: length, ...}, ...}
    :param node_subset: the disjoint subsets of nodes
    :param k: the k shortest paths
    :return:
    """
    # Step 1. Initialization
    neighbor = find_neighbor(network)  # the neighbor set
    v = find_speed(network, neighbor)  # the ripple-spreading speed
    init_node = node_subset[0]
    init_radius = [0]
    init_time = [0]
    init_path = [node_subset[0]]
    temp_path = {}

    # Step 2. The main loop
    for i in range(len(node_subset) - 1):
        destination = node_subset[i + 1]
        init_node, init_time, init_radius, init_path = subRSA(network, neighbor, init_node, destination, init_time, init_radius, init_path, v, k)

    # Step 3. Sort the results
    result = []
    for path in init_path:
        result.append({
            'path': path,
            'length': cal_cost(network, path),
        })
    return result


if __name__ == '__main__':
    # Example 1
    test_network = {
        0: {1: 3, 2: 3},
        1: {0: 3, 2: 3, 3: 5},
        2: {0: 3, 1: 3, 3: 4},
        3: {1: 5, 2: 4},
    }
    subset = [[0], [1], [3]]
    print(main(test_network, subset, 2))

    # Example 2
    # def generate_network():
    #     x = []  # x坐标
    #     y = []  # y坐标
    #     T = [[0], [11, 12, 18, 19], [29, 30, 36, 37], [48]]
    #     connect_list = [0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 20, 21, 27, 28, 34, 35, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
    #     for i in range(7):
    #         for j in range(7):
    #             x.append(i * 10)
    #             y.append(j * 10)
    #     for i in range(1, 49):
    #         x[i] = x[i] + random.uniform(-2, 2)
    #         y[i] = y[i] + random.uniform(-2, 2)
    #     x1 = []
    #     y1 = []
    #     x2 = []
    #     y2 = []
    #     x3 = []
    #     y3 = []
    #     for i in range(49):
    #         if i not in T[1] and i not in T[2]:
    #             x1.append(x[i])
    #             y1.append(y[i])
    #     for i in range(len(T[1])):
    #         x2.append(x[T[1][i]])
    #         y2.append(y[T[1][i]])
    #     for i in range(len(T[2])):
    #         x3.append(x[T[2][i]])
    #         y3.append(y[T[2][i]])
    #     adjacent_matrix = []
    #     for i in range(49):
    #         adjacent_matrix.append([])
    #         for j in range(49):
    #             if i in connect_list and j in connect_list and math.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 15:
    #                 adjacent_matrix[i].append(1)
    #             else:
    #                 adjacent_matrix[i].append(0)
    #     p1 = 0.7
    #     p2 = 0.05
    #     p3 = 0.03
    #     for i in range(49):
    #         for j in range(49):
    #             if (abs(i - j) == 1 or abs(i - j) == 7) and math.sqrt(
    #                     (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 20:  # 横或竖相连
    #                 if random.random() < p1:
    #                     adjacent_matrix[i][j] = 1
    #                     adjacent_matrix[j][i] = 1
    #             if (abs(i - j) == 8 or abs(i - j) == 6) and math.sqrt(
    #                     (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 30:  # 对角线相连
    #                 if random.random() < p2:
    #                     adjacent_matrix[i][j] = 1
    #                     adjacent_matrix[j][i] = 1
    #             if (abs(i - j) == 2 or abs(i - j) == 14) and math.sqrt((x[i] - x[j]) ** 2 + (
    #                     y[i] - y[j]) ** 2) < 30 and i in connect_list and j in connect_list:  # 两横线或两竖线相连
    #                 if random.random() < p3:
    #                     adjacent_matrix[i][j] = 1
    #                     adjacent_matrix[j][i] = 1
    #     for i in range(1, len(T) - 1):
    #         for j in T[i]:
    #             adjacent_matrix[j][j + 7] = 1
    #             adjacent_matrix[j + 7][j] = 1
    #             adjacent_matrix[j][j - 7] = 1
    #             adjacent_matrix[j - 7][j] = 1
    #             adjacent_matrix[j][j + 1] = 1
    #             adjacent_matrix[j + 1][j] = 1
    #             adjacent_matrix[j][j - 1] = 1
    #             adjacent_matrix[j - 1][j] = 1
    #     network = {}
    #     for i in range(49):
    #         network[i] = {}
    #         for j in range(49):
    #             if i != j and adjacent_matrix[i][j] != 0:
    #                 temp_dist = math.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2)
    #                 network[i][j] = temp_dist
    #     return x, y, network
    #
    #
    # def draw_path(x, y, network, subset, path, length, ind):
    #     x1 = []
    #     y1 = []
    #     x2 = []
    #     y2 = []
    #     x3 = []
    #     y3 = []
    #     for i in range(len(x)):
    #         if i in subset[1]:
    #             x2.append(x[i])
    #             y2.append(y[i])
    #         elif i in subset[2]:
    #             x3.append(x[i])
    #             y3.append(y[i])
    #         else:
    #             x1.append(x[i])
    #             y1.append(y[i])
    #     plt.figure(dpi=600)
    #     for i in range(48):
    #         for j in range(i, 49):
    #             if j in network[i].keys():
    #                 temp_x = [x[i], x[j]]
    #                 temp_y = [y[i], y[j]]
    #                 i_index = 49
    #                 j_index = 0
    #                 if i in path and j in path:
    #                     i_index = path.index(i)
    #                     j_index = path.index(j)
    #                 if abs(i_index - j_index) == 1 or (j_index - 1 >= 0 and path[j_index - 1] == i) or (
    #                         i_index + 1 < len(path) and path[i_index + 1] == j):
    #                     plt.plot(temp_x, temp_y, 'black', linewidth=2)
    #                 else:
    #                     plt.plot(temp_x, temp_y, 'springgreen', linewidth=2)
    #     plt.scatter(x1, y1, c='springgreen', s=150)
    #     plt.scatter(x2, y2, c='darkred', alpha=1, s=150)
    #     plt.scatter(x3, y3, c='darkblue', alpha=1, s=150)
    #     plt.title(str(ind) + ', length = ' + str(round(length, 5)))
    #     plt.xticks(())
    #     plt.yticks(())
    #     plt.savefig('D://' + str(ind) + '.png')
    #     plt.show()
    #
    # x, y, test_network = generate_network()
    # subset = [[0], [11, 12, 18, 19], [29, 30, 36, 37], [48]]
    # result = main(test_network, subset, 8)
    # print(result)
    # for i in range(len(result)):
    #     draw_path(x, y, test_network, subset, result[i]['path'], result[i]['length'], i + 1)
