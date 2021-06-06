import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np

NUM_PROCS = 256
TOTAL_TIME = 0.0085

TIME_TICK_PRECISION = 0.000001
TOTAL_TIME_TICKS = int(TOTAL_TIME / TIME_TICK_PRECISION)

with open("heat_equation_result_idle.txt") as file_in:
    lines = file_in.readlines()
    
    matrix = np.ones((NUM_PROCS, TOTAL_TIME_TICKS))
    for line in lines:
        rank, t, it, wait_start, wait_end = line.split();
        start = int(float(wait_start) / TIME_TICK_PRECISION)
        end = int(float(wait_end) / TIME_TICK_PRECISION)
        for i in range(start, end):
            matrix[int(rank)][i] = 0

    x = np.arange(0, TOTAL_TIME_TICKS + 1, 1)
    y = np.arange(0, NUM_PROCS + 1, 1)

    cmap = plt.cm.viridis
    norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)

    im = plt.pcolormesh(x, y, matrix, cmap=cmap, norm=norm, edgecolor="none")
    cbar = plt.colorbar(ticks=[0, 1], shrink=0.2, aspect=8)
    cbar.ax.set_yticklabels(["wait", "active"])
    plt.xlabel("Time (us)")
    plt.ylabel("Rank")
    plt.yticks(range(0, 256 + 1, 16))
    plt.title("Idle period propagation visualization")
    plt.show()
