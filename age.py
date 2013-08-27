import nedclass
import matplotlib.pyplot as plt

nd = nedclass.ned()

zmax = 1500.0
z = 0.001

fp = open('age.out','w')

zplt = []
aplt = []
while z < zmax:
    nd.calcUniverse(z)
    s = '%f\t%f\n' % (z,nd.zage_Gyr)
    fp.write(s)
    zplt.append(z)
    aplt.append(nd.zage_Gyr)
    z = 1.2*z

fp.close()

plt.figure(1)
plt.semilogx(zplt,aplt)
plt.figure(2)
plt.semilogx(aplt,zplt)
plt.figure(3)
plt.loglog(aplt,zplt)
