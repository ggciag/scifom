import numpy as np
import matplotlib.pyplot as plt

topo = np.loadtxt("../ETOPO10.txt")

topo = np.reshape(topo,(1801,3601))

plt.figure(figsize=(20,10))
plt.imshow(topo)
plt.savefig("topo10.png",dpi=300)
plt.close()
