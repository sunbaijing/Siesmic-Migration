import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 6:
    print("\nERROR!, missing params.\nExpected:\n> ./wavefields.py filename x z nt preferred\n")
    exit()
else:
    print("\nDone\n")


data = np.fromfile(sys.argv[1], dtype=np.float32)
data = data.reshape((int(sys.argv[4]),int(sys.argv[2]), int(sys.argv[3])))

plt.matshow(data[int(sys.argv[5]),:,:].T, cmap='seismic')
plt.colorbar()
plt.grid()
plt.title('Wavefield: '+sys.argv[1])
plt.xlabel('X')
plt.ylabel('Depth')
plt.show()
