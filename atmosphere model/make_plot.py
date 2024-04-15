import matplotlib.pyplot as plt

if __name__ == "__main__":
    with open("points.txt", "r") as f:
        f.readline()
        x = f.readline().split(",")
        f.readline()
        f.readline()
        y = f.readline().split(",")
        x.pop()
        y.pop()
        print(x)
        print(y)
