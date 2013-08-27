import math

def calcTCO(MCOmin, z, f_duty=0.1, J=2):
    """Equation 13 from Lidz et al 2011"""
    fcoll = calc_fcoll(MCOmin,z)
    TCO = 2.1*(fcoll/0.1)*(f_duty/0.1)*pow(2.0/J,3.0)*math.sqrt( (1.0+z)/8.0 )
    return TCO
def calc_fcoll(MCOmin, z):
    return 0.1

# 'All' I need to do is implement these...
def calcMoments(MCOmin,k):
    """Equation 16 in Lidz et al"""
    M = 1.0  # This is <M>
    M2 = 1.0/1.0 # This is <M2>
    return M, M2
def calcBias(MCOmin,k):
    """Equation 15 in Lidz et al"""
    b = 1.0/12.0  # This is <b>
    return b
