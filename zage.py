import numpy as np
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/daviddeboer1/Documents/Code/cosmo')
import nedclass as ned
n = ned.ned()
from scipy import interpolate

rest = 1420.40575177
age_of_the_universe = 13.8  #taken from Planck as seen at wikipedia

#redshifts = np.arange(0.001,1000,1)
start = math.log10(0.01)
stop = math.log10(1101.0)
redshifts = np.logspace(start,stop,500)
timeDirection = 'forward'  # vs 'backward'
if timeDirection == 'forward':
    useAge = age_of_the_universe
    dire = -1
    fill_value = 13.
else:
    useAge = 0.0
    dire = -1
    fill_value = 0.0



freq = []
age = []
freq = rest/(1.0+redshifts)
for z in redshifts:
    v = n.calcUniverse(z)
    age.append(useAge - n.DTT_Gyr)
f = interpolate.interp1d(redshifts,age,bounds_error=True,fill_value=fill_value)

#set up x-axis (z,age)
ztics = [1,3,6,15,30,50,100,500]
agetics = f(ztics)
zLabels = []
ageLabels = []
for i,z in enumerate(ztics):
    #print z, agetics[i]
    zLabels.append('%.0f'%(z))
    ageLabels.append('%.2f'%(agetics[i]))
zlim = [1091.0,0.4]
#set up y-axis (freq)
ftics = [30,50,100,150,200,300,500,700]
fLabels = []
for fr in ftics:
    fLabels.append('%.0f'%(fr))

#make plot
plt.figure('FreqOut')
plt.rcParams.update({'font.size':16})
plt.plot(redshifts,freq,linewidth=3)
#eorband stuff
g = interpolate.interp1d(freq[::dire],redshifts[::dire],bounds_error=True,fill_value=0.0)  # the ::-1 just reverses the array since I need increasing array
eorband = [50,100.0,200.0,225]
eorz = g(eorband)
eorages = f(eorz)
eorcolor = ['b','k','k','b']
print '\n\teor band'
print '%6s  %6s  %9s  %7s' % ('MHz ','z  ','Myr  ','Gyr ')
for i,e in enumerate(eorband):
    #z1 = g([eorband[i]])[0]
    print '%6.1f  %6.2f  %9.3f  %7.3f' % (eorband[i],eorz[i],eorages[i]*1000.0,age_of_the_universe-eorages[i])
    clr = eorcolor[i]+'--'
    plt.plot([eorz[i],0.01],[eorband[i],eorband[i]],clr,linewidth=3)

#bottom&right axes
ax1 = plt.subplot(111)
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position('right')
plt.xscale('log')
plt.yscale('log')
plt.xlim(zlim)
plt.ylim([30,1000])
plt.xlabel('Redshift')
plt.ylabel('Frequency [MHz]')
plt.grid()
plt.xticks(ztics,zLabels)
plt.yticks(ftics,fLabels)

#top axis
plt.twiny()
plt.xscale('log')
plt.xlim(zlim)
plt.xticks(ztics,ageLabels)
plt.xlabel('Age [Gyr]')


##plt.xlabel('Age [Gyr]')
##plt.ylabel('Frequency [MHz]')
##plt.semilogx(age,freq)
##
##plt.twiny()
##plt.plot(redshifts,freq)
##plt.xlabel('Redshift')

##plt.rcParams.update({'font.size':16})
###plt.semilogx(age,redshifts)
##ax1 = plt.figure().add_subplot(111)  ##<== this is the magic command in terminal mode to update e.g. xscale
##ax1.plot(age,freq)
##ax1.set_xlabel('Time [Gyr]')
##ax1.set_ylabel('Frequency [MHz]')
##ax1.set_xscale('log')
##ax1.grid()
##ax1.set_xlim([0.05,2])
##ax1.set_ylim([0,300])
##ax2 = ax1.twiny()
###ax2.set_xscale('log')
##ax2.set_xticks(xtics)
##ax2.set_xticklabels(strZlabels)
##ax2.grid(axis='x')
##ax2.set_xlabel('Redshift [z]')
##ax2.plot(redshifts,freq)
###ax2.set_xlim(agelim)

usePlot2 = False
if usePlot2:
    ###PLOT 2:  stuff vs redshift
    option = 2
    f = interpolate.interp1d(age,redshifts,bounds_error=False,fill_value=0.0)
    if option==1:
        ageLabels = [12.6,13.0,13.2,13.4,13.5,13.6]
    elif option==2:
        ageLabels = [12.5,13.5,13.7,13.7]
    elif option==3:
        ageLabels = np.log10(np.logspace(12.5,13.8,6))
    elif option==4:
        ageLabels = [0.1,0.3,0.5,1.0,2.0]
        
    xtics = f(ageLabels)
    for i,x in enumerate(xtics):
        print x,ageLabels[i]
    strAgeLabels = []
    for a in ageLabels:
        strAgeLabels.append(str(a))

    if option == 1:
        zlim = [4,17]
    elif option == 2:
        zlim = [4,27]
    elif option == 3:
        zlim = [4,27]
    elif option == 4:
        zlim = [min(redshifts),max(redshifts)]
    plt.rcParams.update({'font.size':20})
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(redshifts,freq)
    #ax1.set_xlim(zlim)
    #ax1.set_ylim(50,225)
    ax1.set_xlabel('Redshift')
    ax1.set_ylabel('Frequency [MHz]')
    ax1.grid()

    ax2 = ax1.twiny()
    ax2.set_xlabel('Look-back Time [Gyr]')
    ax2.set_xlim(zlim)
    ax2.set_xticks(xtics)
    ax2.set_xticklabels(strAgeLabels)
    ax2.grid(axis='x')

    f = interpolate.interp1d(freq[::dire],redshifts[::dire],bounds_error=False,fill_value=0.0)  # the ::-1 just reverses the array since I need increasing array
    eorband = [50.1,100.0,200.0,224.9]
    eorcolor = ['b','k','k','b']
    for i,e in enumerate(eorband):
        z1 = f([eorband[i]])[0]
        clr = eorcolor[i]+'--'
        ax1.plot([4,z1],[eorband[i],eorband[i]],clr,linewidth=4)

