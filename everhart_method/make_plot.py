import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np

if __name__ == "__main__":
    os.environ["XDG_SESSION_TYPE"] = "xcb"
    with open("points.txt", "r") as f:
        all_file = f.read().split("\n")
        xs = []
        ys = []
        zs = []
        for line in all_file:
            if line != "":
                x, y, z = map(float, line.split(",")[:-1])
                xs.append(x)
                ys.append(y)
                zs.append(z)
    print(len(xs))
    fig = plt.figure()
    ax_3d = fig.add_subplot(projection='3d')
    ax_3d.plot(xs, ys, zs)
    plt.show()
