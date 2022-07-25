import math
import numpy
from scipy import integrate
import transfer
_AP= 0.322
_aP= 0.75
_qP= 0.3
def tophat_window(k,R):
    x= k*R
    return 3.*(numpy.sin(x)-x*numpy.cos(x))/x**3.

def S(R,zend,vbc,ncost=21,full=False,short=True,window='tophat',
      dm=None,dmks=None):
    """
    NAME:
       S
    PURPOSE:
       calculate S=\sigma(R)^2, fluctuations within a sphere of radius R
    INPUT:
       R - radius (Mpc)
       zend - redshift
       vbc - velocity
       ncost= number of cost to use
       full=, short=
       window= window function to use
       dm= matter transfer function to use (if provided)
       dmks= ks for the dm
    OUTPUT:
       S(R)
    HISTORY:
       2012-03-08 - Written - Bovy (IAS)
    """
    #Calculate matter power spectrum
    if dm is None:
        dmks= transfer.ks
        nks= len(dmks)
        dm= numpy.zeros(nks)
        for ii in range(nks):
            #BOVY: KOSHER? Don't think it matters much bc this is all high mass
            #BOVY: ACTUALLY NOT USED, BC DM IS PRECALCULATED
            if ks[ii] < 1.: full, short= False, False
            else: full, short= True, False
            dm[ii]= transfer.dmAvg(zend,ks[ii],vbc,ncost=21,full=full,short=short)
    if window.lower() == 'tophat':
        return integrate.trapz((dmks/0.002)**transfer._NSMINUSONE*dm**2.*dmks**-1.*transfer._DELTA2ZETA*tophat_window(dmks,R),dmks)

def dndm(M,zend,vbc,ncost=21,full=False,short=True,dc=1.67,thisS=None,
         ps=False):
    #Sheth-Tormen / Press Schechter mass function
    # Use as: dndm= massfunction.dndm(20.,0.,ncost=1,short=False)
    #   calculates this for all ms= transfer.kfToMF(transfer.ks)
    #Jo: Figure out units! Units are /Mpc^3
    #Calculate S
    if thisS is None:
        thisS= numpy.zeros(len(M))
        for ii in range(len(M)):
            thisS[ii]= S(transfer.mToR(M[ii]),
                         zend,vbc,ncost=21,full=False,short=True)
    dSdm= (thisS-numpy.roll(thisS,1))/(M-numpy.roll(M,1))
    dSdm[0]= dSdm[1] #Hack that does not matter (these masses are like 10^26!)
    dSdm= numpy.fabs(dSdm)
    if ps:
        return transfer._RHOO/M*dSdm*fps(thisS,dc=dc)
    else:
        return transfer._RHOO/M*dSdm*fst(thisS,dc=dc)

def fst(S,dc=1.67):
    #Sheth and Tormen (1999) function
    nu= dc/numpy.sqrt(S)
    return _AP*nu/S*numpy.sqrt(_aP/2./math.pi)*(1.+1./(_aP*nu**2.)**_qP)\
        *numpy.exp(-_aP*nu**2./2.)
def fps(S,dc=1.67):
    #Original Press-Schechter function
    #Re-definition, does not change globally
    _aP= 1.
    _qP= 0.
    _AP= 0.5
    nu= dc/numpy.sqrt(S)
    return _AP*nu/S*numpy.sqrt(_aP/2./math.pi)*(1.+1./(_aP*nu**2.)**_qP)\
        *numpy.exp(-_aP*nu**2./2.)

def dmstardm(mh,zend,vbc,ncost=21,fs=0.01,thismf=None,rlss=None,
             **kwargs):
    if thismf is None or rlss is None:
        thismf,rlss= transfer.mf(zend,vbc,ncost=ncost,returnrlss=True,**kwargs)
    thisfgas= transfer.fgas(mh,thismf,rlss)
    return fs*thisfgas+fs*3.*transfer.fb0(rlss)*(2.**(transfer._ALPHAFGAS/3.)-1.)*(thismf/mh)**transfer._ALPHAFGAS*(1.+(2.**(transfer._ALPHAFGAS/3.)-1.)*(thismf/mh)**transfer._ALPHAFGAS)**(-3./transfer._ALPHAFGAS-1.)

def dmvdmstar(mstar):
    return 2.5/math.log(10.)/mstar
