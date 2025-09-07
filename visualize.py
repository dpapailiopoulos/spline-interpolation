import sys
import matplotlib.pyplot as plt

try:
    path = sys.argv[1]
except IndexError:
    path = "./build/cubic_spline_points.txt"

f = open(path, 'r')

X = []
Y = []
for row in f:
    row = row.strip().split()
    X.append(float(row[0]))
    Y.append(float(row[1]))
fig, ax = plt.subplots()
fig.set_size_inches(13,9)
ax.plot(X, Y, color='orange')
ax.axis('equal')

plt.show()