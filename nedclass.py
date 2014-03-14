from math import *

class ned:
  def __init__(self,z=6.0, H0=68.0, WM=0.3086, WV=0.6914):
    """Cosmology calculator ala Ned Wright (www.astro.ucla.edu/~wright)
      input values = redshift, Ho, Omega_m, Omega_vac
      ouput values = age at z, distance in Mpc, kpc/arcsec, apparent to abs mag conversion"""

    print 'H0, WM, WV taken from Planck values as reported at wikipedia'
    self.z = z
    self.H0 = H0
    self.h = self.H0/100.0
    h = self.h
    self.WM = WM
    if WV==None:
      self.WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda
    else:
      self.WV = WV
    self.WR = 4.165E-5/(h*h)    # Omega(radiation) includes 3 massless neutrino species, T0 = 2.72528    
    self.WK = 1-self.WM-self.WR-self.WV        # Omega curvaturve = 1-Omega(total)
    self.c = 299792.458         # velocity of light in km/sec
    self.Tyr = 977.8            # coefficent for converting 1/H into Gyr
    self.DTT = 0.5              # time from z to now in units of 1/H0
    self.DTT_Gyr = 0.0          # value of DTT in Gyr
    self.age = 0.5              # age of Universe in units of 1/H0
    self.age_Gyr = 0.0          # value of age in Gyr
    self.zage = 0.1             # age of Universe at redshift z in units of 1/H0
    self.zage_Gyr = 0.0         # value of zage in Gyr
    self.DCMR = 0.0             # comoving radial distance in units of c/H0
    self.DCMR_Mpc = 0.0 
    self.DCMR_Gyr = 0.0
    self.DA = 0.0               # angular size distance
    self.DA_Mpc = 0.0
    self.DA_Gyr = 0.0
    self.kpc_DA = 0.0
    self.DL = 0.0               # luminosity distance
    self.DL_Mpc = 0.0
    self.DL_Gyr = 0.0           # DL in units of billions of light years
    self.V_Gpc = 0.0
    self.az = 1.0/(1+1.0*z)     # 1/(1+z(object))

    print "z = "+str(self.z)
    print "H0 = "+str(self.H0)
    print "WM = "+str(self.WM)
    print "WV = "+str(self.WV)
    print "WR = "+str(self.WR)

  def calcUniverse(self,z,returnValue = 'DA_Mpc', verbose=False,n=1000):
    """Calculates cosmological values"""
    H0 = self.H0
    h = self.h
    WR = self.WR
    WM = self.WM
    WV = self.WV
    WK = self.WK
    c = self.c
    Tyr = self.Tyr
    az = 1.0/(1.0+z)

    self.age = 0.0
    for i in range(n):
      a = az*(i+0.5)/n
      adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
      self.age = self.age + 1./adot
    self.zage = az*self.age/n
    self.zage_Gyr = (Tyr/H0)*self.zage

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    self.DTT = 0.0
    self.DCMR = 0.0
    for i in range(n):
      a = az+(1-az)*(i+0.5)/n
      adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
      self.DTT = self.DTT + 1./adot
      self.DCMR = self.DCMR + 1./(a*adot)

    self.DTT = (1.-az)*self.DTT/n
    self.DCMR = (1.-az)*self.DCMR/n
    self.age = self.DTT+self.zage
    self.age_Gyr = self.age*(Tyr/H0)
    self.DTT_Gyr = (Tyr/H0)*self.DTT
    self.DCMR_Gyr = (Tyr/H0)*self.DCMR
    self.DCMR_Mpc = (c/H0)*self.DCMR
    self.Ez = sqrt(self.WM*pow(1.0+z,3.0) + self.WV)
    self.Ynu = (3000.0/self.Ez/self.h)*pow(1.0+z,2 )  # this is Y*restFreq

    DCMR = self.DCMR
    # tangential comoving distance
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
      if WK > 0:
        ratio =  0.5*(exp(x)-exp(-x))/x 
      else:
        ratio = sin(x)/x
    else:
      y = x*x
      if WK < 0: y = -y
      ratio = 1. + y/6. + y*y/120.

    self.DCMT = ratio*DCMR
    self.DA = az*self.DCMT
    self.DA_Mpc = (c/H0)*self.DA
    self.kpc_DA = self.DA_Mpc/206.264806
    self.DA_Gyr = (Tyr/H0)*self.DA
    self.DL = self.DA/(az*az)
    self.DL_Mpc = (c/H0)*self.DL
    self.DL_Gyr = (Tyr/H0)*self.DL
    self.X = (1.0+z)*self.DA_Mpc

    # comoving volume computation
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
      if WK > 0:
        ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
      else:
        ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
    else:
      y = x*x
      if WK < 0: y = -y
      ratio = 1. + y/5. + (2./105.)*y*y

    self.VCM = ratio*DCMR*DCMR*DCMR/3.
    self.V_Gpc = 4.*pi*((0.001*c/H0)**3)*self.VCM

    if verbose:
      print 'For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ',
      print '%1.2f' % WV + ', z = ' + '%1.3f' % z
      print 'It is now ' + '%1.1f' % self.age_Gyr + ' Gyr since the Big Bang.'
      print 'The age at redshift z was ' + '%1.1f' % self.zage_Gyr + ' Gyr.'
      print 'The light travel time was ' + '%1.1f' % self.DTT_Gyr + ' Gyr.'
      print 'The comoving radial distance, which goes into Hubbles law, is',
      print '%1.1f' % self.DCMR_Mpc + ' Mpc or ' + '%1.1f' % self.DCMR_Gyr + ' Gly.'
      print 'The comoving volume within redshift z is ' + '%1.1f' % self.V_Gpc + ' Gpc^3.'
      print 'The angular size distance D_A is ' + '%1.1f' % self.DA_Mpc + ' Mpc or',
      print '%1.1f' % self.DA_Gyr + ' Gly.'
      print 'comoving DA is '+str(self.DA_Mpc*(1.0+z))
      print 'This gives a scale of ' + '%.2f' % self.kpc_DA + ' kpc/".'
      print 'The luminosity distance D_L is ' + '%1.1f' % self.DL_Mpc + ' Mpc or ' + '%1.1f' % self.DL_Gyr + ' Gly.'
      print 'The distance modulus, m-M, is '+'%1.2f' % (5*log10(self.DL_Mpc*1e6)-5)
##    else:
##      print '%1.2f' % self.zage_Gyr,
##      print '%1.2f' % self.DCMR_Mpc,
##      print '%1.2f' % self.kpc_DA,
##      print '%1.2f' % (5*log10(self.DL_Mpc*1e6)-5)

    ##print 'z = ',z

    if (returnValue == 'DA_Mpc'):
      val = self.DA_Mpc
    elif (returnValue == 'Ez'):
      val = self.Ez
    else:
      val = self.DA_Mpc
      
    return val

