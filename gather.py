
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 4:
    print("\nERROR!, missing params.\nExpected:\n> ./gather.py filename x time\n")
    exit()
else:
    print("\nDone\n")


data = np.fromfile(sys.argv[1], dtype=np.float32)
data = data.reshape((int(sys.argv[2]),int(sys.argv[3])))

plt.matshow(data.T, cmap='seismic')
plt.colorbar()
plt.grid()
plt.title('Gather: '+sys.argv[1])
plt.xlabel('X')
plt.ylabel('Time')
plt.show()
