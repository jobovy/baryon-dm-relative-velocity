import math
import numpy
from scipy import integrate, optimize, interpolate
_OMEGAB= 0.0456
_OMEGAC= 0.227
_OMEGAL= 0.728
_ZEQ= 3232.
_H0NOC= 70.4
_DELTA2ZETA= 2.44*10.**-9.
_OMEGAM= _OMEGAB+_OMEGAC
_FB= _OMEGAB/_OMEGAM
_FC= _OMEGAC/_OMEGAM
_SQRTOMEGAM= math.sqrt(_OMEGAM)
_ZRECOMB= 1020.
_TCMB= 2.726
_A1= 1./119.
_A2= 1./115.
_GAMMACS= 5./3.
_GAMMAKTB0OVERMUMH= _GAMMACS*8.617332*_TCMB/1.22/931.*10.**-11.
_KTB0OVERMUMH= 8.617332*_TCMB/1.22/931.*10.**-11.
_KOVERMUMH= 8.617332/1.22/931.*10.**-11.
_C= 300000.
_H0= _H0NOC/_C
_NSMINUSONE= -0.037
_RHOO= 3./8./math.pi*_H0**2./4.302/10.**-9.*_C**2.*_OMEGAM
_XEINIT= 6.5*10**-2 #1.2*10.**-5*math.sqrt(_OMEGAM)/_OMEGAB/(_HONOC/100.)
_PARSECKM= 3.08568025*10.**13.
_tGAMMAINV= 8.55*10.**-13./(365.24*24*60*60.)/_C*10.**6.*_PARSECKM
_ALPHABCOEFF=2.753*10.**-14./_C*10.**6.*_PARSECKM #cm^3/Mpc
_TWOTHI= 2.*157807.
_YHE=0.079
_NHZCOEFF= 8.6*10.**-6.*_OMEGAB*(_H0NOC/100.)**2.
#_NHZCOEFF= 3.*_H0NOC**2./8./math.pi/(4.302)*10.**-6.*1.989/1.66/(3.08)**3.*_OMEGAB #This is the value that Cora used
_ALPHAFGAS= 0.7
#Load the initial conditions
def readInit():
    """Read the initial conditions from CAMB"""
    Deltas= numpy.loadtxt('../tex/Delta_z%i.dat' % int(_ZRECOMB))
    ks= Deltas[:,0]
    Deltac= Deltas[:,1]
    Deltab= Deltas[:,2]
    Deltaph= Deltas[:,3]
    dDeltab= numpy.loadtxt('../tex/dDeltabdz_z%i.dat' % int(_ZRECOMB))[:,1]
    dDeltac= numpy.loadtxt('../tex/dDeltacdz_z%i.dat' % int(_ZRECOMB))[:,1]
    dDeltaph= numpy.loadtxt('../tex/dDeltardz_z%i.dat' % int(_ZRECOMB))[:,1]
    return (ks,Deltac,Deltab,dDeltac,dDeltab,Deltaph,dDeltaph)
ks,Deltac,Deltab,dDeltac,dDeltab, Deltaph, dDeltaph= readInit()
ksOne= ks[(ks > 1.)]
#Read CAMB xe
#xe_recfast= numpy.loadtxt('../tex/xe_recfast.dat')
_NOREIONIZATION= False
if _NOREIONIZATION:
    xe_recfast= numpy.loadtxt('../tex/xe_without_reion.dat')
else:
    xe_recfast= numpy.loadtxt('../tex/xe_with_reion.dat')
indx= (xe_recfast[:,0] < _ZRECOMB+10.)
xe_recfast= xe_recfast[indx,:]
if _NOREIONIZATION:
    logxe_rec= interpolate.UnivariateSpline(xe_recfast[:,0][::-1],
                                            numpy.log(xe_recfast[:,1])[::-1])
else:
    logxe_rec= interpolate.interp1d(xe_recfast[:,0][::-1],
                                    numpy.log(xe_recfast[:,1])[::-1],
                                    kind='linear',bounds_error=False,
                                    fill_value=numpy.log(xe_recfast[-1,1]))
#xe_jo= numpy.loadtxt('../output/xe_z1020.dat',delimiter='|')

#Read CAMB cs2
cs_camb= numpy.loadtxt('../output/z_cs_tb.dat')
indx= (cs_camb[:,0] >= 9.)
cs_z= cs_camb[indx,0]
cs2_camb= cs_camb[indx,1]
cs2_camb_interp= interpolate.UnivariateSpline(numpy.log(cs_z),numpy.log(cs2_camb))
def logxe_jo(z):
    if isinstance(z,numpy.ndarray):
        indx= (numpy.floor((_ZRECOMB-z)/(xe_jo[0,0]-xe_jo[1,0]))).astype('int')
    else:
        indx= int(numpy.floor((_ZRECOMB-z)/(xe_jo[0,0]-xe_jo[1,0])))
    return numpy.log(xe_jo[indx,1])

def matterHz(z):
    """H(z) for matter-domination, although not anymore!"""
    return _H0*numpy.sqrt(_OMEGAM*(1.+z)**3.*(1.+(1.+z)/(1.+_ZEQ))+_OMEGAL)
#return _H0*_SQRTOMEGAM*(1.+z)**1.5*(1.+(1.+z)/(1.+_ZEQ))**0.5

#dlogks= numpy.roll(numpy.log(ks),-1)-numpy.log(ks)
#dlogks[-1]= dlogks[-2]
#sbc= matterHz(_ZRECOMB)*(1.+_ZRECOMB)**0.\
#    *math.sqrt(_DELTA2ZETA*numpy.sum((dDeltab-dDeltac)**2./ks**2.*dlogks))
sbc= matterHz(_ZRECOMB)\
    *math.sqrt(_DELTA2ZETA*integrate.trapz((ks/0.002)**_NSMINUSONE*(dDeltab-dDeltac)**2./ks**3.,ks))

def tbato(a):
    return 1./a/(1.+a/_A1/(1.+(_A2/a)**1.5))

def dlntbatodlna(a):
    return -(4.*a**4 + 2.* _A1*_A2**3 + 4.*a**3.*_A1*(_A2/a)**1.5+a**3.*(2.*_A1+7.*_A2*(_A2/a)**.5))/(2.*a*(a + _A2*(_A2/a)**.5)*(a**2.+a*_A1+_A1*_A2*(_A2/a)**.5))

def cs2(z):
    #return _GAMMAKTB0OVERMUMH*tbato(1./(1.+z))
    return _KTB0OVERMUMH*tbato(1./(1.+z))*(1.-1./3.*dlntbatodlna(1./(1.+z)))

def alphaB(T):
    """'case B' recombination coeff in cm^3/Mpc"""
    lamb= _TWOTHI/T
    return _ALPHABCOEFF*lamb**1.5/(1.+(lamb/2.740)**0.407)**2.242

def zrecomb():
    return _ZRECOMB

def kInit():
    return ks

def fb0(rlss):
    return _OMEGAB/_OMEGAM*(1+3.2*rlss)

def fgas(m,mf,rlss):
    return fb0(rlss)*(1.+(2.**(_ALPHAFGAS/3.)-1.)*(mf/m)**_ALPHAFGAS)**(-3./_ALPHAFGAS)

def fgasnew(m,mf,rlss,beta=0.97,gamma=2.15):
    if beta is None:
        beta= _ALPHAFGAS
    if gamma is None:
        gamma= beta/3.
    return fb0(rlss)*(1.+(2.**gamma-1.)*(mf/m)**beta)**(-1./gamma)

def mstar(mh,zend,vbc,ncost=21,fs=0.01,thismf=None,rlss=None,new=False,
          **kwargs):
    """Calculate Mstar = Mhalo = fs * fgas(zend,vbc,mh,...), use this to convert halo mass into stellar mass"""
    #First calculate mf for this vbc and zend
    if thismf is None or rlss is None:
        thismf,rlss= mf(zend,vbc,ncost=ncost,returnrlss=True,new=new,
                        **kwargs)
    if new:
        thisfgas= fgasnew(mh,thismf,rlss)
    else:
        thisfgas= fgas(mh,thismf,rlss)
    return mh*fs*thisfgas

def mv(mh,zend,vbc,ncost=21,fs=0.01,mvpersolarmass=6.7,new=False,
       thismstar=None,**kwargs):
    if thismstar is None:
        thismstar= mstar(mh,zend,vbc,ncost=ncost,fs=fs,new=new,**kwargs)
    return mvpersolarmass-2.5*numpy.log10(thismstar)

def kfToMF(kff):
    return 4./3.*math.pi*_RHOO*(math.pi/kff)**3.
def mTok(m):
    return (3.*m/4./math.pi/_RHOO)**(-1./3.)*math.pi
def mToR(m):
    return (3.*m/4./math.pi/_RHOO)**(1./3.)

def vcirc(mhalo,z):
    """Circular velocity from Koposov et al. (2009), in km/s"""
    return (4.302*10.**-9.*mhalo/rvir(mhalo,z))**0.5
def rvir(mhalo,z):
    """Virial radius from Koposov et al. (2009), in Mpc"""
    return (3.*mhalo/(4.*math.pi*_rvirDeltac(z)*_RHOO*(1.+z)**3.))**(1./3.)
def _rvirDeltac(z):
    x= -(1.-_OMEGAM)/(_OMEGAM*(1.+z)**3.+1.-_OMEGAM)
    return (18.*math.pi**2.+82.*x-39.*x**2.)/(1.+x)

def mf(zend,vbc,ncost=21,alt=True,fitdbdm=False,fixrlss=True,full=False,
       short=True,returnrlss=False,new=False):
    kff= kfAvg(zend,vbc,ncost=ncost,alt=alt,fitdbdm=fitdbdm,fixrlss=fixrlss,
               full=full,short=short,returnrlss=returnrlss)
    if new:
        kff[0]= kff[0]/numpy.sqrt(1.+numpy.fabs(vbc)/sbc)
    if returnrlss: return (kfToMF(kff[0]),kff[1])
    else: return kfToMF(kff)

def mfAvg(zend,vbc,ncost=21,fixrlss=True,alt=True,full=False,
          short=True,returnrlss=False):
    costs= numpy.linspace(-1.,1.,ncost)
    out= 0.
    if returnrlss: outrlss= 0.
    for cost in costs:
        kff= kf(zend,vbc,cost,full=full,fixrlss=fixrlss,alt=alt,
                short=short,returnrlss=returnrlss)
        if returnrlss:
            out+= kfToMF(kff[0])
            outrlss+= kff[1]
        else:
            out+= kfToMF(kff)
    if returnrlss: return (out/ncost,outrlss/ncost)
    else: return out/ncost

def kf(zend,vbc,cost,full=False,alt=True,fixrlss=True,short=True,
       returnrlss=False):
    #Calculate dbdtotAvg for each k > 1 Mpc^-1
    dbdt= numpy.zeros(len(ksOne))
    for ii in range(len(ksOne)):
        dbdt[ii]= dbdtot(zend,ksOne[ii],vbc,cost=cost,full=full,short=short)
    #Now fit for kf
    #indx= (dbdt > 0.3)
    #kf= fitKF(ksOne[indx],dbdt[indx],alt=alt,fixrlss=fixrlss)
    kf= fitKF(ksOne,dbdt,alt=alt,fixrlss=fixrlss,returnrlss=returnrlss)
    return kf

def kfAvg(zend,vbc,ncost=21,alt=True,fitdbdm=False,fixrlss=True,full=False,
          short=True,returnrlss=False):
    #Calculate dbdtotAvg for each k > 1 Mpc^-1
    dbdt= numpy.zeros(len(ksOne))
    for ii in range(len(ksOne)):
        if fitdbdm:
            dbdt[ii]= dbdmAvg(zend,ksOne[ii],vbc,ncost=ncost,full=full,
                              short=short)
        else:
            dbdt[ii]= dbdtotAvg(zend,ksOne[ii],vbc,ncost=ncost,full=full,
                                short=short)
    #Now fit for kf
    #indx= (dbdt > 0.3)
    #kf= fitKF(ksOne[indx],dbdt[indx],alt=alt,fixrlss=fixrlss)
    kf= fitKF(ksOne,dbdt,alt=alt,fixrlss=fixrlss,returnrlss=returnrlss)
    return kf
#return (kf,ksOne,dbdt)

def fitKF(k,dbdt,alt=True,fixrlss=True,returnrlss=False):
    """Fit 1-k^2/kF^2+rlss"""
    if isinstance(fixrlss,bool) and fixrlss:
        smallkIndx= (k > 1.)*(k < 10.)
        smallk= k[smallkIndx]
        smallkdbdt= dbdt[smallkIndx]
        rlss= numpy.mean(smallkdbdt)-1.
        print rlss
    elif isinstance(fixrlss,float):
        rlss= fixrlss
        fixrlss= True
    else:
        rlss= None
    if alt:
        if fixrlss:
            init_params= [numpy.log(400.),numpy.log(10.)]
            _func= _fitKFFuncAltFixedrlss
        else:
            init_params= [numpy.log(400.),-0.05,10.]
            _func= _fitKFFuncAlt
    else:
        if fixrlss:
            init_params= [numpy.log(400.)]
            _func= _fitKFFuncFixedrlss
        else:
            init_params= [numpy.log(400.),-0.05]
            _func= _fitKFFunc
    params= optimize.fmin_powell(_func,init_params,(k,dbdt,rlss))
    print params
    if not fixrlss: rlss= params[1]
    try:
        if returnrlss: return (numpy.exp(params[0]),rlss)
        else: return numpy.exp(params[0])
    except IndexError:
        if returnrlss: return (numpy.exp(params),rlss)
        else: return numpy.exp(params)

def _fitKFFunc(p,k,dbdt,rlss):
    return numpy.sum((dbdt-1.-p[1]+k**2./numpy.exp(p[0])**2.)**2.)

def _fitKFFuncAlt(p,k,dbdt,rlss):
    pred= (1.+p[1])*numpy.power(1.+1/p[2]*k**2./numpy.exp(p[0])**2./(1.+p[1]),-p[2])
    return numpy.sum((dbdt-pred)**2.)

def _fitKFFuncFixedrlss(p,k,dbdt,rlss):
    return numpy.sum((dbdt-1.-rlss+k**2./numpy.exp(p[0])**2.)**2.)

def _fitKFFuncAltFixedrlss(p,k,dbdt,rlss):
    pred= (1.+rlss)*numpy.power(1.+1/numpy.exp(p[1])*k**2./numpy.exp(p[0])**2./(1.+rlss),-numpy.exp(p[1]))
    return numpy.sum((dbdt-pred)**2.)

def dmAvg(zend,k,vbc,ncost=21,full=False,short=True):
    #Calculate dm for all of these cost
    costs= numpy.linspace(-1.,1.,ncost)
    out= 0.
    for cost in costs:
        out+= dm(zend,k,vbc,cost,full=full,short=short)
    return out/ncost

def dbdtotAvg(zend,k,vbc,ncost=21,full=False,short=True):
    #Calculate dbdtot for all of these cost
    costs= numpy.linspace(-1.,1.,ncost)
    out= 0.
    for cost in costs:
        out+= dbdtot(zend,k,vbc,cost,full=full,short=short)
    return out/ncost

def dbdmAvg(zend,k,vbc,ncost=21,full=False,short=True):
    #Calculate dbdm for all of these cost
    costs= numpy.linspace(-1.,1.,ncost)
    out= 0.
    for cost in costs:
        out+= dbdm(zend,k,vbc,cost,full=full,short=short)
    return out/ncost

def dbdtot(zend,k,vbc,cost,full=False,short=True):
    #Solve transfer eqns
    if short:
        yout= solveTransferEqsShort(zend,k,vbc,cost)
        return numpy.sqrt(yout[2]**2.+yout[7]**2.)/\
            (_FB*numpy.sqrt(yout[2]**2.+yout[7]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[5]**2.))
    elif full:
        yout= solveTransferEqsFull(zend,k,vbc,cost)
        return numpy.sqrt(yout[2]**2.+yout[9]**2.)/\
            (_FB*numpy.sqrt(yout[2]**2.+yout[9]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[7]**2.))
    else:
        yout= solveTransferEqsComplex(zend,k,vbc,cost)
        return numpy.sqrt(yout[2]**2.+yout[6]**2.)/\
            (_FB*numpy.sqrt(yout[2]**2.+yout[6]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[4]**2.))

def dm(zend,k,vbc,cost,full=False,short=True):
    #Solve transfer eqns
    if short:
        yout= solveTransferEqsShort(zend,k,vbc,cost)
        return (_FB*numpy.sqrt(yout[2]**2.+yout[7]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[5]**2.))
    elif full:
        yout= solveTransferEqsFull(zend,k,vbc,cost)
        return (_FB*numpy.sqrt(yout[2]**2.+yout[9]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[7]**2.))
    else:
        yout= solveTransferEqsComplex(zend,k,vbc,cost)
        return (_FB*numpy.sqrt(yout[2]**2.+yout[6]**2.)+_FC*numpy.sqrt(yout[0]**2.+yout[4]**2.))

def dTdb(zend,k,vbc,cost):
    #Solve transfer eqns
    yout= solveTransferEqsShort(zend,k,vbc,cost)
    return numpy.sqrt(yout[4]**2.+yout[9]**2.)/numpy.sqrt(yout[2]**2.+yout[7]**2.)

def dbdm(zend,k,vbc,cost,full=False,short=True):
    #Solve transfer eqns
    if short:
        yout= solveTransferEqsShort(zend,k,vbc,cost)
    elif full:
        yout= solveTransferEqsFull(zend,k,vbc,cost)
    else:
        yout= solveTransferEqs(zend,k,vbc,cost)
    return yout[2]/(yout[0]+yout[2])

def initdeltaT(deltag,ddeltag,thetab):
    return deltag/4.-((-matterHz(_ZRECOMB)*(1.+_ZRECOMB)*ddeltag/4.)
                      +2./3.*thetab)/_tGAMMAINV/_XEINIT/(1.+_ZRECOMB)**4.

def solveTransferEqs(zend,k,vbc,cost):
    #Find the initial conditions for this k
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx],
                     Deltab[indx],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx]]).reshape(4)
    return odeintTransferEqs(yo,zend,vbc,k,cost)

def solveTransferEqsComplexAvg(zend,k,vbc,ncost=11):
    costs= numpy.linspace(-1.,1.,ncost)
    out= numpy.zeros((8,ncost))
    for ii in range(ncost):
        out[:,ii]= solveTransferEqsComplex(zend,k,vbc,costs[ii])
    return numpy.sum(out,axis=1)/ncost

def solveTransferEqsComplex(zend,k,vbc,cost):
    #Find the initial conditions for this k
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx][0],
                     Deltab[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx][0],
                     0.,
                     0.,
                     0.,
                     0.]).reshape(8)
    return odeintTransferEqsComplex(yo,zend,vbc,k,cost)

def xez(zend):
    #Find the initial conditions for this k
    k= ks[0] #just pick a k
    vbc, cost= 0., 0.
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx],
                     Deltab[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx],
                     Deltaph[indx][0]/4.,
                     _XEINIT,
                     _TCMB*(1.+_ZRECOMB)]).reshape(7)
    t= numpy.linspace(zrecomb(),zend,10001)
    out= integrate.odeint(transferEqsFull,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[:,5]

def Tz(zend):
    #Find the initial conditions for this k
    k= ks[0] #just pick a k
    vbc, cost= 0., 0.
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx],
                     Deltab[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx],
#                     Deltaph[indx][0]/4.,
                     initdeltaT(Deltaph[indx][0],dDeltaph[indx]/4.,
                                matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx]),
                     _XEINIT,
                     _TCMB*(1.+_ZRECOMB)]).reshape(7)
    t= numpy.linspace(zrecomb(),zend,1001)
    out= integrate.odeint(transferEqsFull,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[:,6]


def solveTransferEqsFull(zend,k,vbc,cost):
    #Find the initial conditions for this k
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx],
                     Deltab[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx],
                     initdeltaT(Deltaph[indx][0],dDeltaph[indx]/4.,
                                matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx]),
#                     Deltaph[indx][0]/4.,
                     _XEINIT,
                     _TCMB*(1.+_ZRECOMB),
                     0.,
                     0.,
                     0.,
                     0.,
                     0.]).reshape(12)
    return odeintTransferEqsFull(yo,zend,vbc,k,cost)

def solveTransferEqsShort(zend,k,vbc,cost):
    #Find the initial conditions for this k
    indx= (ks == k)
    if numpy.sum(indx) == 0.:
        raise IOError("Requested k not available")
    yo= numpy.array([Deltac[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltac[indx],
                     Deltab[indx][0],
                     matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx],
                     #Deltaph[indx][0]/4.,
                     initdeltaT(Deltaph[indx][0],dDeltaph[indx]/4.,
                                matterHz(_ZRECOMB)*(1.+_ZRECOMB)*dDeltab[indx]),
                     0.,
                     0.,
                     0.,
                     0.,
                     0.]).reshape(10)
                     
    return odeintTransferEqsShort(yo,zend,vbc,k,cost)

def odeintTransferEqs(yo,zend,vbc,k,cost):
    #t= numpy.array([zrecomb(),zend])
    t= numpy.linspace(zrecomb(),zend,1001)
    out= integrate.odeint(transferEqs,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[-1,:]

def odeintTransferEqsComplex(yo,zend,vbc,k,cost):
    #t= numpy.array([zrecomb(),zend])
    t= numpy.linspace(zrecomb(),zend,1001)
    out= integrate.odeint(transferEqsComplex,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[-1,:]

def transferEqs(y,z,vbc,k,cost):
    """Transfer equations RHS"""
    thisH= matterHz(z)
    a= (1.+_ZRECOMB)/(1.+z)
    thisvbc= vbc/a
    #Sound speed
    thiscs2= cs2(z) #_GAMMAKTB0OVERMUMH*tbato(1./(1.+z))
    #cs2= numpy.exp(cs2_camb_interp(numpy.log(z)))
    #print thisvbc/math.sqrt(cs2)
    a= 1./(1.+z)
    return -numpy.array([cost*thisvbc*k/thisH*y[0]-y[1]/thisH*a,
                         -cost*thisvbc*k/thisH*y[1]
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                         -2.*a*y[1],
                         -y[3]/thisH*a,
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                         -2.*a*y[3]+thiscs2*k**2./a/thisH*y[2]])

def transferEqsComplex(y,z,vbc,k,cost):
    """Transfer equations RHS"""
    thisH= matterHz(z)
    a= (1.+_ZRECOMB)/(1.+z)
    thisvbc= vbc/a
    #Sound speed
    thiscs2= cs2(z) #_GAMMAKTB0OVERMUMH*tbato(1./(1.+z))
    #cs2= numpy.exp(cs2_camb_interp(numpy.log(z)))
    #print thisvbc/math.sqrt(cs2)
    a= 1./(1.+z)
    return -numpy.array([-cost*thisvbc*k/thisH*y[4]-y[1]/thisH*a,
                          -cost*thisvbc*k/thisH*y[5]
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                          -2.*a*y[1],
                          -y[3]/thisH*a,
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                          -2.*a*y[3]+thiscs2*k**2./a/thisH*y[2],
                          #imaginary part
                          cost*thisvbc*k/thisH*y[0]-y[5]/thisH*a,
                          cost*thisvbc*k/thisH*y[1]
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[4]+_FB*y[6])
                          -2.*a*y[5],
                          -y[7]/thisH*a,
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[4]+_FB*y[6])
                          -2.*a*y[7]+thiscs2*k**2./a/thisH*y[6]])
#    return -numpy.array([-y[1]/thisH*a,
#                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
#                          -2.*a*y[1],
#                          cost*thisvbc*k/thisH*y[6]-y[3]/thisH*a,
#                          cost*thisvbc*k/thisH*y[7]
#                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
#                          -2.*a*y[3]+thiscs2*k**2./a/thisH*y[2],
#                          #imaginary part
#                          -y[5]/thisH*a,
#                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[4]+_FB*y[6])
#                          -2.*a*y[5],
#                          -cost*thisvbc*k/thisH*y[2]-y[7]/thisH*a,
#                          -cost*thisvbc*k/thisH*y[3]
#                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[4]+_FB*y[6])
#                          -2.*a*y[7]+thiscs2*k**2./a/thisH*y[6]])

def odeintTransferEqsShort(yo,zend,vbc,k,cost):
    #t= numpy.array([zrecomb(),zend])
    t= numpy.linspace(zrecomb(),zend,1001)
    out= integrate.odeint(transferEqsShort,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[-1,:]

def transferEqsShort(y,z,vbc,k,cost):
    """Transfer equations RHS, including spatially varying sound speed, but using pre-computed xe and Tb"""
    thisH= matterHz(z)
    a= (1.+_ZRECOMB)/(1.+z)
    thisvbc= vbc/a
    #print thisvbc/math.sqrt(cs2)
    a= 1./(1.+z)
    tcmb= _TCMB/a
    nhz= _NHZCOEFF/a**3.
    #if z >= 11.:
    #    xe= numpy.exp(logxe_rec(z))
    #else:
    #    xe= 1.
    xe= numpy.exp(logxe_rec(z))
    tb= _TCMB*tbato(1./(1.+z))
    #if z < 20.: tb*= 1.1
    return -numpy.array([-cost*thisvbc*k/thisH*y[5]-y[1]/thisH*a,
                          -cost*thisvbc*k/thisH*y[6]
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                          -2.*a*y[1],
                          -y[3]/thisH*a,
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                          -2.*a*y[3]
                          +_KOVERMUMH*tb*(y[2]+y[4])*k**2./a/thisH,
                          -2./3.*y[3]*a/thisH
                          -xe*_tGAMMAINV/a**3./thisH*tcmb/tb*y[4],
                          #Imaginary part
                          cost*thisvbc*k/thisH*y[0]-y[6]/thisH*a,
                          cost*thisvbc*k/thisH*y[1]
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[5]+_FB*y[7])
                          -2.*a*y[6],
                          -y[8]/thisH*a,
                          -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[5]+_FB*y[7])
                          -2.*a*y[8]
                          +_KOVERMUMH*tb*(y[7]+y[9])*k**2./a/thisH,
                          -2./3.*y[8]*a/thisH
                         -xe*_tGAMMAINV/a**3./thisH*tcmb/tb*y[9]])

def odeintTransferEqsFull(yo,zend,vbc,k,cost):
    #t= numpy.array([zrecomb(),zend])
    t= numpy.linspace(zrecomb(),zend,1001)
    out= integrate.odeint(transferEqsFull,yo,t,args=(vbc,k,cost),
                          rtol=10.**-8.)
    return out[-1,:]

def transferEqsFull(y,z,vbc,k,cost):
    """Transfer equations RHS, including spatially varying sound speed"""
    thisH= matterHz(z)
    a= (1.+_ZRECOMB)/(1.+z)
    thisvbc= vbc/a
    #print thisvbc/math.sqrt(cs2)
    a= 1./(1.+z)
    tcmb= _TCMB/a
    nhz= _NHZCOEFF/a**3.
    return -numpy.array([-cost*thisvbc*k/thisH*y[7]-y[1]/thisH*a,
                          -cost*thisvbc*k/thisH*y[8]
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                         -2.*a*y[1],
                         -y[3]/thisH*a,
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[0]+_FB*y[2])
                         -2.*a*y[3]
                          +_KOVERMUMH*y[6]*(y[2]+y[4])*k**2./a/thisH,
                         -2./3.*y[3]*a/thisH
                         -y[5]*_tGAMMAINV/a**3./thisH*tcmb/y[6]*y[4],
                         -alphaB(y[6])*y[5]**2.*(1.+_YHE)*nhz*a/thisH,
                         -2.*a*y[6]+y[5]*_tGAMMAINV*(tcmb-y[6])/a**3/thisH,
                          #Imaginary part, only first 5 equations
                         cost*thisvbc*k/thisH*y[0]-y[8]/thisH*a,
                         -cost*thisvbc*k/thisH*y[1]
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[7]+_FB*y[9])
                         -2.*a*y[8],
                         -y[10]/thisH*a,
                         -1.5/thisH*_H0**2.*_OMEGAM/a**2.*(_FC*y[7]+_FB*y[9])
                         -2.*a*y[10]
                          +_KOVERMUMH*y[6]*(y[9]+y[11])*k**2./a/thisH,
                          -2./3.*y[10]*a/thisH
                          -y[5]*_tGAMMAINV/a**3./thisH*tcmb/y[6]*y[11]])
