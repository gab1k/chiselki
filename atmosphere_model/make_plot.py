import matplotlib.pyplot as plt

if __name__ == "__main__":
    with open("points_logs.txt", "r") as f:
        f.readline()
        x = f.readline().split(",")
        f.readline()
        f.readline()
        y = f.readline().split(",")
        x.pop()
        y.pop()
        x = list(map(float, x))
        y = list(map(float, y))
        print(x)
        print(y)
        print(len(x), len(y))
        plt.scatter(x, y)
        plt.title("траектория полета снаряда")
        plt.xlabel("Длина (м)")
        plt.ylabel("Высота (м)")
        plt.show()
