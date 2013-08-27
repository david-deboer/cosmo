import subprocess
import math
import matplotlib.pyplot as plt
import cosmoFunc as cfnc

class powerSpectrum:
    def __init__(self, z=3.0, h=0.7, n = 0.95, k=0.1,MCOmin=1.0E8,f_duty=0.1,J=2):
        self.z = z
        self.h = h
        self.n = n
        if type(k) != list:
            self.k = [k]
        else:
            self.k = k
        self.MCOmin = MCOmin
        self.f_duty = f_duty
        self.J = J

    def linDelta2(self, k=None,z=None, h=None, n=None):
        """P lin"""
        if z == None:
            z = self.z
        else:
            self.z = z
        if h == None:
            h = self.h
        else:
            self.h = h
        if n == None:
            n = self.n
        else:
            self.n = n
        if k == None:
            k = self.k
        else:
            if type(k)!=list:
                k = [k]
            self.k = k
        TF = self.transferFunction(z,h,k)
        d_H = self.delta_H(n)
        Ho = 100.0*h
        c = 3.0E5  # in km/sec
        self.D2lin = []
        for i,ki in enumerate(k):
            self.D2lin.append( d_H*d_H * pow(c*ki/Ho,3.0+n) * TF[i]*TF[i] )
        return self.D2lin

    def cluster(self,k=None,z=None):
        if z==None:
            z = self.z
        if k==None:
            k = self.k
        TCO = cfnc.calcTCO(self.MCOmin, z, self.f_duty, self.J)
        b = []
        for ki in k:
            b.append(cfnc.calcBias(self.MCOmin,ki))
        D2lin = self.linDelta2(k=k,z=z)
        D2cluster = []
        for i,v in enumerate(D2lin):
            D2cluster.append(TCO*TCO*b[i]*b[i]*v)
        self.D2cluster = D2cluster
        return D2cluster

    def poisson(self,k=None):
        if k==None:
            k = self.k
        TCO = cfnc.calcTCO(self.MCOmin, self.z, self.f_duty, self.J)
        M=[]
        M2=[]
        for ki in k:
            mi, m2i = cfnc.calcMoments(self.MCOmin,ki)
            M.append(mi)
            M2.append(m2i)
        D2poisson = []
        for i in range(len(M)):
            D2poisson.append(pow(k[i],3.0)/(2.0*math.pi*math.pi)*TCO*TCO*M2[i]/self.f_duty/(M[i]*M[i]))
        self.D2poisson = D2poisson
        return D2poisson

    def powerSpectrum(self,k=None,z=None):
        if k==None:
            k=self.k
        if z==None:
            z=self.z
        C = self.cluster(k,z)
        P = self.poisson(k)
        D2 = []
        for i in range(len(C)):
            D2.append(C[i] + P[i])
        self.D2 = D2
        return D2

    def plotP(self,z=None,h=None,n=None,kStart=1.0E-4,dlogk=0.05,kStop=10.0,save2self=True):
        if z==None:
            z = self.z
        if h==None:
            h = self.h
        if n == None:
            n = self.n
        k = genk(kStart,dlogk,kStop)
        D2 = self.powerSpectrum(k,z)

        plt.loglog(k,D2,label='total')
        plt.loglog(k,self.D2cluster,label='cluster')
        plt.loglog(k,self.D2poisson,label='poisson')
        plt.xlabel('$k$ [Mpc$^{-1}$]')
        plt.ylabel('$\Delta^2(k,z)$')
        plt.title('Power Spectrum z=%.1f, h=%.2f' % (z,h))
        plt.legend(loc='lower right')
        if (save2self):
            self.k = k
            self.D2 = D2


    def plotPlin(self,z=None,h=None,n=None,kStart=1.0E-4,dlogk=0.05,kStop=10.0,save2self=True):
        if z==None:
            z = self.z
        if h==None:
            h = self.h
        if n == None:
            n = self.n
        k = genk(kStart,dlogk,kStop)
        D2lin = self.linDelta2(k,z,h,n)

        plt.loglog(k,D2lin)
        plt.xlabel('$k$ [Mpc$^{-1}$]')
        plt.ylabel('$\Delta^2_{lin}(k,z)$')
        plt.title('Linear Power Spectrum z=%.1f, h=%.2f' % (z,h))
        if (save2self):
            self.k = k
            self.D2lin = D2lin

    def transferFunction(self, z=None, h=None, k=None, tmpFile='trans.out'):
        """Calls c function trans (from transferFunc.c) and returns value(s)"""
        if z == None:
            z = self.z
        else:
            self.z = 6.0
        if h == None:
            h = self.h
        else:
            self.h = h
        n = self.n
        if k == None:
            k = self.k
        else:
            if type(k)!=list:
                k = [k]
            self.k = k

        self.TF = []
        for ki in k:
            comm = "./trans %f %f %f %f > %s" % (z,h,n,ki,tmpFile)
            subprocess.call(comm, shell=True)
            fp = open(tmpFile,'r')
            line = fp.readline()
            if line[0]=='T':
                data = float(line.split()[1])
            else:
                data = 0.0
            self.TF.append(data)
            #print '==>  ', ki, data
            fp.close()

            
        return self.TF

    def plotTransferFunction(self,z=None,h=None,kStart=1.0E-4,dlogk=0.05,kStop=10.0,save2self=True):
        if z==None:
            z = self.z
        if h==None:
            h = self.h
        n = self.n
        k = genk(kStart,dlogk,kStop)
        TF = self.transferFunction(z,h,k)

        plt.loglog(k,TF)
        plt.xlabel('k [Mpc$^{-1}$]')
        plt.ylabel('T(k)')
        plt.title('Transfer function z=%.1f, h=%.2f' % (z,h))
        if (save2self):
            self.k = k
            self.TF = TF

    def delta_H(self,n=None,tmpFile='delta.out'):
        if n==None:
            n = self.n
        else:
            self.n = n
        z = self.z
        h = self.h
        ki = 1.0
        comm = "./trans %f %f %f %f > %s" % (z,h,n,ki,tmpFile)
        subprocess.call(comm, shell=True)
        fp = open(tmpFile,'r')
        line = fp.readline()
        line = fp.readline()
        if line[0]=='D':
            self.delta = float(line.split()[1])
        else:
            self.delta = 0.0
        return self.delta

def genk(kStart=1.0E-4,dlogk=0.02,kStop=10.0):
    k = []
    kStart = math.log10(kStart)
    kStop = math.log10(kStop)
    klog = kStart
    while klog <= kStop:
        k.append(pow(10.0,klog))
        klog+=dlogk
    return k
    
