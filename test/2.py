import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("1.txt")[:20]
print(data)

fig = plt.figure()
ax = fig.add_subplot(projection="polar")
c = ax.scatter(data[:, 1], data[:, 0], alpha=0.75)
plt.show()
