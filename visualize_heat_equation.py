import matplotlib.pyplot as plt
import numpy as np

with open("heat_equation_result.txt") as file_in:
    matrix = []

    first_line = file_in.readline().split();
    matrix.append(first_line)
    dim = len(first_line)
    for _ in range(dim - 1):
        matrix.append(file_in.readline().split())

    fig, ax = plt.subplots()
    r = np.arange(0, dim + 1, 1)
    matrix = np.array(matrix).astype(np.float)
    ax.pcolormesh(r, r, matrix, cmap="hot")
    ax.axis("equal")
    plt.show()
