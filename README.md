### The Ripple-Spreading Algorithm for the k Shortest Path Problem

##### The shortest path tour problem (SPTP) aims to find the shortest path that traverses multiple disjoint node subsets in a given order. The *k*-SPTP aims to find the k shortest paths for the SPTP.

----

| Variables     | Meaning                                                      |
| ------------- | ------------------------------------------------------------ |
| network       | Dictionary, {node1: {node2: length, node3: length, ...}, ...} |
| node_subset   | List, [[subset1], [subset2], ...]                            |
| source        | List, the source nodes of this subproblem of SPTP            |
| destination   | List, the destination nodes of this subproblem of SPTP       |
| init_time     | List, the initial time that should generate initial ripples at source nodes |
| init_radius   | List, the initial radius of initial ripples at source nodes  |
| init_path     | List, the initial path for each initial ripple               |
| k             | The k shortest paths                                         |
| nn            | The number of nodes                                          |
| neighbor      | Dictionary, {node1: [the neighbor nodes of node1], ...}      |
| v             | The ripple-spreading speed (i.e., the minimum length of arcs) |
| t             | The simulated time index                                     |
| nr            | The number of ripples - 1                                    |
| epicenter_set | List, the epicenter node of the ith ripple is epicenter_set[i] |
| path_set      | List, the path of the ith ripple from the source node to node i is path_set[i] |
| radius_set    | List, the radius of the ith ripple is radius_set[i]          |
| active_set    | List, active_set contains all active ripples                 |
| Omega         | Dictionary, Omega[n] = i denotes that ripple i is generated at node n |

----

#### Example 1

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/k-SPTP%20example.png)

```python
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
```

##### Output

```python
[
    {'path': [0, 1, 3], 'length': 8}, 
    {'path': [0, 1, 2, 3], 'length': 10}
]
```

---------

#### Example 2

##### This example aims to find the 8 shortest paths for the *k*-SPTP on a randomly generated network with 49 nodes. 

```python
if __name__ == '__main__':
    # Example 2
    def generate_network():
        x = []  # x坐标
        y = []  # y坐标
        T = [[0], [11, 12, 18, 19], [29, 30, 36, 37], [48]]
        connect_list = [0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 20, 21, 27, 28, 34, 35, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
        for i in range(7):
            for j in range(7):
                x.append(i * 10)
                y.append(j * 10)
        for i in range(1, 49):
            x[i] = x[i] + random.uniform(-2, 2)
            y[i] = y[i] + random.uniform(-2, 2)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        x3 = []
        y3 = []
        for i in range(49):
            if i not in T[1] and i not in T[2]:
                x1.append(x[i])
                y1.append(y[i])
        for i in range(len(T[1])):
            x2.append(x[T[1][i]])
            y2.append(y[T[1][i]])
        for i in range(len(T[2])):
            x3.append(x[T[2][i]])
            y3.append(y[T[2][i]])
        adjacent_matrix = []
        for i in range(49):
            adjacent_matrix.append([])
            for j in range(49):
                if i in connect_list and j in connect_list and math.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 15:
                    adjacent_matrix[i].append(1)
                else:
                    adjacent_matrix[i].append(0)
        p1 = 0.7
        p2 = 0.05
        p3 = 0.03
        for i in range(49):
            for j in range(49):
                if (abs(i - j) == 1 or abs(i - j) == 7) and math.sqrt(
                        (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 20:  # 横或竖相连
                    if random.random() < p1:
                        adjacent_matrix[i][j] = 1
                        adjacent_matrix[j][i] = 1
                if (abs(i - j) == 8 or abs(i - j) == 6) and math.sqrt(
                        (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) < 30:  # 对角线相连
                    if random.random() < p2:
                        adjacent_matrix[i][j] = 1
                        adjacent_matrix[j][i] = 1
                if (abs(i - j) == 2 or abs(i - j) == 14) and math.sqrt((x[i] - x[j]) ** 2 + (
                        y[i] - y[j]) ** 2) < 30 and i in connect_list and j in connect_list:  # 两横线或两竖线相连
                    if random.random() < p3:
                        adjacent_matrix[i][j] = 1
                        adjacent_matrix[j][i] = 1
        for i in range(1, len(T) - 1):
            for j in T[i]:
                adjacent_matrix[j][j + 7] = 1
                adjacent_matrix[j + 7][j] = 1
                adjacent_matrix[j][j - 7] = 1
                adjacent_matrix[j - 7][j] = 1
                adjacent_matrix[j][j + 1] = 1
                adjacent_matrix[j + 1][j] = 1
                adjacent_matrix[j][j - 1] = 1
                adjacent_matrix[j - 1][j] = 1
        network = {}
        for i in range(49):
            network[i] = {}
            for j in range(49):
                if i != j and adjacent_matrix[i][j] != 0:
                    temp_dist = math.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2)
                    network[i][j] = temp_dist
        return x, y, network


    def draw_path(x, y, network, subset, path, length, ind):
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        x3 = []
        y3 = []
        for i in range(len(x)):
            if i in subset[1]:
                x2.append(x[i])
                y2.append(y[i])
            elif i in subset[2]:
                x3.append(x[i])
                y3.append(y[i])
            else:
                x1.append(x[i])
                y1.append(y[i])
        plt.figure(dpi=600)
        for i in range(48):
            for j in range(i, 49):
                if j in network[i].keys():
                    temp_x = [x[i], x[j]]
                    temp_y = [y[i], y[j]]
                    i_index = 49
                    j_index = 0
                    if i in path and j in path:
                        i_index = path.index(i)
                        j_index = path.index(j)
                    if abs(i_index - j_index) == 1 or (j_index - 1 >= 0 and path[j_index - 1] == i) or (
                            i_index + 1 < len(path) and path[i_index + 1] == j):
                        plt.plot(temp_x, temp_y, 'black', linewidth=2)
                    else:
                        plt.plot(temp_x, temp_y, 'springgreen', linewidth=2)
        plt.scatter(x1, y1, c='springgreen', s=150)
        plt.scatter(x2, y2, c='darkred', alpha=1, s=150)
        plt.scatter(x3, y3, c='darkblue', alpha=1, s=150)
        plt.title(str(ind) + ', length = ' + str(round(length, 5)))
        plt.xticks(())
        plt.yticks(())
        plt.savefig('D://' + str(ind) + '.png')
        plt.show()

    x, y, test_network = generate_network()
    subset = [[0], [11, 12, 18, 19], [29, 30, 36, 37], [48]]
    result = main(test_network, subset, 8)
    print(result)
    for i in range(len(result)):
        draw_path(x, y, test_network, subset, result[i]['path'], result[i]['length'], i + 1)
```

##### Output

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/1.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/2.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/3.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/4.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/5.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/6.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/7.png)

![](https://github.com/Xavier-MaYiMing/The-ripple-spreading-algorithm-for-the-k-shortest-path-tour-problem/blob/main/8.png)

```python
[
    {'path': [0, 1, 2, 10, 17, 18, 17, 24, 23, 30, 31, 32, 33, 34, 41, 48], 'length': 150.1495608090849}, 
    {'path': [0, 1, 2, 10, 17, 18, 17, 24, 23, 30, 37, 44, 45, 46, 47, 48], 'length': 150.63572706262775}, 
    {'path': [0, 1, 2, 10, 17, 18, 25, 32, 31, 30, 31, 32, 33, 34, 41, 48], 'length': 151.2016145482568}, 
    {'path': [0, 1, 8, 9, 17, 18, 17, 24, 23, 30, 31, 32, 33, 34, 41, 48], 'length': 151.27653543201578}, 
    {'path': [0, 1, 2, 10, 17, 18, 17, 24, 23, 30, 31, 32, 39, 40, 47, 48], 'length': 151.29119265478346}, 
    {'path': [0, 1, 2, 10, 11, 10, 17, 24, 23, 30, 31, 32, 33, 34, 41, 48], 'length': 151.44287618906438}, 
    {'path': [0, 1, 2, 10, 17, 18, 25, 32, 31, 30, 37, 44, 45, 46, 47, 48], 'length': 151.68778080179965}, 
    {'path': [0, 1, 8, 9, 17, 18, 17, 24, 23, 30, 37, 44, 45, 46, 47, 48], 'length': 151.76270168555862}
]
```

