import nedclass
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
sys.path.append('/Users/ddeboer/Documents/Code/Lib')
import fitlib

n = nedclass.ned()
restFreq = 230.0   # GHz
BW = 1.0           # GHz
nCh = 1024
zi = 7

zList = [3.0, 6.0]
nurList = [115.0, 230.0]

###Make fitCosmo1
for restFreq in nurList:
    z = []  # redshift
    X = []  # X Mpc/rad
    Xb = [] # scaled version
    Y = []  # Y Mpc/GHz
    DA = [] # angular scale term
    E = []  # Hubble parameter evolution
    apX = []    # approximated X (post-fit numbers put back in)
    apY = []    # approximated Y ( " )
    z2fit = []  # range of z's to fit for approximation
    X2fit = []  #     "    X's   "
    Y2fit = []  #     "    Y's   "

    # Loop over red-shift
    for i in range(100):
        zi = i*11./99.
        z.append(zi)
        DAi = n.calcUniverse(zi)
        DA.append(DAi)
        Xi = DAi*(1.+zi)
        obsFreq = restFreq/(1.+zi)
        obsWavelength = 0.3/obsFreq
        print 'z=%.2f:  rest = %.1f GHz, obs = %.2f GHz (%.4f m)' % (zi,restFreq,obsFreq,obsWavelength)
        Xbi = Xi/((180.0/math.pi)*1.0)  #convert to degrees
        X.append(Xi)
        Xb.append(Xbi)
        apX.append(4670.*pow(1.+zi,0.312))
        E.append(n.Ez)
        Yi = ( n.Ynu )/restFreq
        Y.append(Yi)
        apY.append(64.*pow(1.+zi,0.56))
        if zi > 6.0 and zi < 9.0:
            z2fit.append(zi)
            X2fit.append(Xi)
            Y2fit.append(Yi)

    if restFreq == 115.0:
        lineColor = 'b'
    else:
        lineColor = 'g'
    lineType = '-'
    lineChar = lineColor+lineType
    if restFreq == nurList[0]:  #(only do once)
        plt.figure(1)
        #plt.subplot(211)
        plt.plot(z,X,lineChar)
        plt.xlabel('Z')
        plt.ylabel('X [Mpc/rad]')#  and  Y [Mpc/GHz]')
        plt.plot(z,apX)

        ##fp = open('X.dat','w')
        ##for i,zi in enumerate(z):
        ##    s = '%f\t%f\n' %(zi,X[i])
        ##    fp.write(s)
        ##fp.close()

        z2fit = np.array(z2fit)
        X2fit = np.array(X2fit)
        fitlib.fit(z2fit, X2fit,showPlot=True)

    lineType = '--'
    lineChar = lineColor+lineType
    #plt.figure(2)
    #plt.subplot(212)
    #plt.plot(z,Y,lineChar)
    #plt.xlabel('Z')
    #plt.ylabel('Y [Mpc/GHz]')
    ##plt.plot(z,apY)
    ##Y2fit = np.array(Y2fit)
    ##fitlib.fit(z2fit,Y2fit,showPlot=False)


###Make "fitCosmo10"
###Perpendicular
for zi in zList:
    for restFreq in nurList:
        DAi = n.calcUniverse(zi)
        Xi = DAi*(1.+zi)
        obsFreq = restFreq/(1.+zi)
        obsWavelength = 0.3/obsFreq  # m
        nbstep = 100
        bmin = 0.5
        bmax = 10.0
        bstep = (bmax-bmin)/nbstep
        Mpcperp = []
        b = []
        kperp = []

        for i in range(nbstep):
            bi = bmin + i*bstep
            ki = 2.0*math.pi*bi/(obsWavelength*Xi)
            Mpci = 2.0*math.pi/ki
            b.append(bi)
            Mpcperp.append(Mpci)
            kperp.append(ki)

        plt.figure(10)
        if restFreq == 115.0:
            lineColor = 'b'
        else:
            lineColor = 'g'
        if zi == 3.0:
            lineType = '--'
        else:
            lineType = '-'
        lineChar = lineColor+lineType

        plt.subplot(221)
        plt.semilogy(b,kperp,lineChar)
        plt.xlabel('Baseline [m]')
        plt.ylabel('k [Mpc$^{-1}$]')
        plt.grid()

        #plt.figure(11)
        plt.subplot(223)
        plt.plot(b,Mpcperp,lineChar)
        plt.xlabel('Baseline [m]')
        plt.ylabel('Size [Mpc]')
        plt.grid()

        ###Parallel
        DAi = n.calcUniverse(zi)
        Yi = n.Ynu/restFreq
        Mpcpar = []
        kpar = []
        bw = []

        for i in range(nCh):
            bi = BW/(i+1)
            bw.append(bi)
            ki = 2.0*math.pi/(bi*Yi)
            kpar.append(ki)
            Mpci = 2.0*math.pi/ki
            Mpcpar.append(Mpci)

        #plt.figure(20)
        plt.subplot(222)
        plt.semilogy(bw,kpar,lineChar)
        plt.xlabel('Band [GHz]')
        plt.ylabel('k [Mpc$^{-1}$]')
        plt.grid()

        #plt.figure(21)
        plt.subplot(224)
        plt.plot(bw,Mpcpar,lineChar)
        plt.xlabel('Band [GHz]')
        plt.ylabel('Size [Mpc]')
        plt.grid()


restFreq = 230.0#/2.
zi = 6.42
#zi = 2.71
obsFreq = restFreq/(1.0+zi)
obsWavelength = 0.3/obsFreq
BW = 31.2/1000.0
BW = 2.
DAi = n.calcUniverse(zi)
Xi = DAi*(1.+zi)
fov = obsWavelength/1.2
spatial = Xi*fov
Yi = n.Ynu/restFreq
ki = 2.0*math.pi/(BW*Yi)
Mpci = 2.0*math.pi/ki
print '----'
print 'X, spatial = ',Xi, spatial
print 'nYnu, DAi, zi: ',n.Ynu,DAi,zi
print '::::  %f  %f kpar (BW=%.3f MHz) = %f\tMpc_par = %f' % (restFreq,restFreq/(1.+zi),BW*1000.0,ki,Mpci)
restFreq = 230.0#/2.
zi = 6.42
#zi = 2.71
obsFreq = restFreq/(1.0+zi)
obsWavelength = 0.3/obsFreq
BW = 31.2/1000.0
BW = 2.
DAi = n.calcUniverse(zi)
Xi = DAi*(1.+zi)
fov = obsWavelength/4.8
spatial = Xi*fov
Yi = n.Ynu/restFreq
ki = 2.0*math.pi/(BW*Yi)
Mpci = 2.0*math.pi/ki
print '----'
print 'X, spatial = ',Xi, spatial
print 'nYnu, DAi, zi: ',n.Ynu,DAi,zi
print '::::  %f  %f kpar (BW=%.3f MHz) = %f\tMpc_par = %f' % (restFreq,restFreq/(1.+zi),BW*1000.0,ki,Mpci)



paperRestFreq = 1.42
zi = 6
Yi = n.Ynu/paperRestFreq
paperBW = 0.008
ki = 2.0*math.pi/(paperBW*Yi)
Mpci = 2.0*math.pi/ki
print 'PAPER:  kpar (BW=%.3f MHz) = %f\tMpc_par = %f' % (paperBW*1000.0,ki,Mpci)
