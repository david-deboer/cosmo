import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/daviddeboer1/Documents/Code/cosmo')
import nedclass as ned
n = ned.ned()
from scipy import interpolate

rest = 1420.40575177

redshifts = np.arange(3,50,.2)

freq = []
age = []
for z in redshifts:
    freq.append(rest/(1.0+z))
    v = n.calcUniverse(z)
    age.append(n.DTT_Gyr)

f = interpolate.interp1d(age,redshifts)
ageLabels = [12.6,12.8,13.0,13.2,13.4,13.5,13.6]
xtics = f(ageLabels)
for i,x in enumerate(xtics):
    print x,ageLabels[i]
strAgeLabels = []
for a in ageLabels:
    strAgeLabels.append(str(a))

zlim = [4,17]
plt.rcParams.update({'font.size':20})
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(redshifts,freq)
ax1.set_xlim(zlim)
ax1.set_ylim(50,225)
ax1.set_xlabel('Redshift')
ax1.set_ylabel('Frequency [MHz]')
ax1.grid()

ax2 = ax1.twiny()
ax2.set_xlabel('Look-back Time [Gyr]')
ax2.set_xlim(zlim)
ax2.set_xticks(xtics)
ax2.set_xticklabels(strAgeLabels)
ax2.grid(axis='x')

f = interpolate.interp1d(freq[::-1],redshifts[::-1])  # the ::-1 just reverses the array since I need increasing array
eorband = [100.0,200.0]
z1 = f([eorband[0]])[0]
z2 = f([eorband[1]])[0]
ax1.plot([4,z1],[eorband[0],eorband[0]],'k--')
ax1.plot([4,z2],[eorband[1],eorband[1]],'k--')

