import numpy as np
import matplotlib.pyplot as plt
E=np.zeros(200*200)
E=E.reshape(200,200)
for i in range(200):
   E[i,int(i*i/800)+160:199]=3100
   E[i,120:int(i*i/800)+160]=2700
E[:,90:120]=2200
E[:,20:90]=1900
for i in range(200):
    E[i,20:int(i/4)]=1200
E[:,0:20]=800
plt.matshow(E.T)
plt.colorbar()
plt.show()
E=E.astype(np.float32)
E.tofile('synthmodel2D.bin')
