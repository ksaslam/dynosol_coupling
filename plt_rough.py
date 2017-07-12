import numpy as np
import matplotlib.pyplot as plt
import seistools
x = np.linspace(0., 80, 5)
y= np.ones(5) * 20.0 + seistools.rough.generate_profile(5 , 80 , 10**-2 , 20 , h = 1 , seed= 200)

# plt.plot(x,y, 'o')
# plt.show()
np.savetxt('x.out', x)
np.savetxt('y.out', y)
