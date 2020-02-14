
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 5:
    print("\nERROR!, missing params.\nExpected:\n> ./ill.py filename z x clip\n")
    exit()
else:
    print("\nDone.\n")


data = np.fromfile(sys.argv[1], dtype=np.float32)
data = data.reshape((int(sys.argv[3]),int(sys.argv[2])))

plt.matshow(data.T, vmax=np.percentile(data, float(sys.argv[4])), cmap='Greys')
plt.colorbar()
plt.grid()
plt.title('Illumination: '+sys.argv[1])
plt.xlabel('X')
plt.ylabel('Depth')
plt.show()
