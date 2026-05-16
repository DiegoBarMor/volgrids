import numpy as np
from collections import deque

benzene = np.array([ # 1SWI
    [-10.248,  -9.558,  31.480],
    [-10.010,  -8.648,  32.511],
    [ -8.779,  -7.983,  32.540],
    [ -7.827,  -8.144,  31.532],
    [ -8.124,  -8.997,  30.470],
    [ -9.315,  -9.722,  30.458],
])

ros = np.array([ # 1F1T
    [111.283, -75.053, -81.421],
    [111.875, -76.278, -81.938],
    [111.580, -76.503, -83.284],
    [112.210, -77.766, -83.761],
    [111.944, -77.958, -84.994],
    [112.592, -79.217, -85.361],
    [112.404, -79.594, -86.738],
    [111.666, -78.803, -87.653],
    [111.039, -77.567, -87.223],
    [111.202, -77.177, -85.844],
    [110.505, -75.888, -85.477],
    [110.813, -75.697, -84.118],
    [110.220, -74.472, -83.577],
    [110.509, -74.216, -82.201],
    [109.762, -75.149, -86.378],
    [110.121, -73.882, -86.894],
    [109.392, -73.130, -87.817],
    [108.135, -73.649, -88.335],
    [107.799, -74.912, -87.791],
    [108.496, -75.698, -86.863],
    [111.487, -74.688, -80.098],
    [112.301, -75.489, -79.188],
    [110.950, -73.485, -79.427],
    [113.013, -80.819, -87.060],
    [113.777, -81.783, -86.237],
    [112.885, -81.300, -88.481],
])

spi = np.array([ # 1D8F
    [ 3.918,  48.681,  53.575],
    [ 3.376,  47.785,  52.598],
    [ 5.348,  48.751,  53.527],
    [ 3.284,  48.486,  55.244],
    [ 3.949,  49.164,  56.343],
    [ 3.579,  48.828,  57.713],
    [ 2.556,  47.826,  57.967],
    [ 1.894,  47.168,  56.853],
    [ 2.248,  47.486,  55.465],
    [ 2.133,  47.478,  59.190],
    [ 2.835,  47.524,  60.324],
    [ 3.641,  50.224,  52.920],
    [ 2.315,  50.311,  52.291],
    [ 2.298,  51.565,  51.409],
    [ 2.849,  52.749,  52.095],
    [ 4.213,  52.608,  52.633],
    [ 4.322,  51.387,  53.520],
    [ 2.259,  53.933,  51.975],
    [ 2.699,  54.923,  52.619],
    [ 1.271,  54.034,  51.207],
    [ 0.805,  55.232,  50.647],
    [-0.387,  55.003,  49.771],
    [-0.991,  53.690,  49.839],
    [-2.135,  53.405,  48.991],
    [-2.598,  54.480,  48.115],
    [-1.954,  55.795,  48.081],
    [-0.802,  56.096,  48.923],
    [ 1.140,  50.302,  53.282],
    [ 0.280,  49.289,  53.214],
    [-0.859,  49.482,  54.031],
    [ 0.954,  51.193,  54.129],
])

MAX_BOND_LENGTH = 2.1

def dfs(graph: np.ndarray, idx_start: int) -> tuple[int]:
    paths  : deque[set] = deque()
    queue  : deque[int] = deque()
    parents: deque[int] = deque()

    paths.append(set())
    queue.append(idx_start)
    parents.append(-1)

    idxs = np.arange(len(graph))
    while queue:
        node = int(queue.popleft())
        path = paths.popleft()
        parent = parents.popleft()

        neighs = [int(i) for i in idxs[graph[node]] if i != parent]

        if len(path) and node == idx_start:
            return tuple(sorted(path))

        if node in path: continue

        path.add(node)

        queue.extend(neighs)
        paths.extend([path.copy() for _ in neighs])
        parents.extend([node for _ in neighs])


def find_aromatics(coords):
    mat0 = coords.reshape(1,-1,3)
    mat1 = coords.reshape(-1,1,3)
    dist = np.linalg.norm(mat1 - mat0, axis = -1)

    bonds = (0 < dist) & (dist < MAX_BOND_LENGTH)

    cycles = set(
        dfs(bonds, i) for i in range(len(coords))
    ) - {None}

    for cycle in cycles:
        if len(cycle) <= 3: continue

        normals = np.zeros((len(cycle), 3))
        for idx,i in enumerate(cycle):
            j = cycle[(idx+1) % len(cycle)]
            k = cycle[(idx+2) % len(cycle)]
            n = np.cross(coords[i] - coords[j], coords[k] - coords[j])
            n /= np.linalg.norm(n)
            normals[idx] = n

        dev = np.linalg.norm(
            np.std(normals, axis = 0)
        )
        print(cycle, "DEVIATION:", dev)


print("benzene"); find_aromatics(benzene)
print("ros"); find_aromatics(ros)
print("spi"); find_aromatics(spi)
