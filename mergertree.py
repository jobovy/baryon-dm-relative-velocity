#
# current usage:
# 
# out= mergertree.mergertree(dz=10.**-5.,zend=11.,dzsub=0.01)
# zs, out, zacc, newhalos, newhalosprerei, zcmain, subzc= out
#
import sys
import math
import numpy
from scipy import integrate
import cosmolopy.perturbation as cp
from transfer import _OMEGAB, _OMEGAM, _OMEGAL, _H0NOC, _NSMINUSONE
import massfunction
import transfer
_ERASESTR= "                                                                                "
_NMSINTEGRATE= 1001
_NCOLLAPSEZS= 1001
cosmo= {'N_nu': 0,
        'Y_He': 0.24,
        'baryonic_effects': False,
        'h': _H0NOC/100.,
        'n': 1-_NSMINUSONE,
        'omega_M_0': _OMEGAM,
        'omega_b_0': _OMEGAB,
        'omega_lambda_0': _OMEGAL,
        'omega_n_0': 0.0,
        'sigma_8': 0.809,
        't_0': 13.75,
        'tau': 0.087,
        'w': -1.0,
        'z_reion': 11.}

def mergertree(minit=10.**12.,zend=11.,mres=10.**5.,dz=10.**-5.,zstart=0.,
               dzsub=10.**-2.,justcollapsez=False,
               dm=None,dmks=None,dontevolvesub=False,
               _nosubstructure=False,_fixdz=False,_zs=None,
               _logmsgtmres=None,_logmsltmres=None,
               _siggtmres=None,_sigltmres=None,
               _dsigdmgtmres=None,_dsigdmltmres=None):
    """
    NAME:
       mergertree
    PURPOSE:
       generate a merger tree, using the Cole et al. (2000) prescription and the Eisenstein & Hu transfer function, unless dm is given
    INPUT:
       minit= (10^12) mass of the main z=0 halo / Msolar
       zend= (11) highest redshift
       mres= (10^5) simulate halos down to this mass / Msolar
       dz= (10.**-5.) initial redshift resolution of the tree
       dzsub= (10.**-3.) redshift resolution for following subhalos
       justcollapsez= (False) just calculate the collapse redshift (useful for testing)
       dm=, dmks= matter transfer function and the ks it was calculated for
       dontevolvesub= (False) If true, don't evolve subhalos
    INTERNAL:
       _nosubstructure= do not merge substructure (for evolving subhalos)
       _fixdz= fix the redshift grid to _zs
       _zs= fixed redshift grid
       +_logmsgtmres etc.
    OUTPUT:
        (zs,array [nz,zacc,newhalo,newhaloprerei) 
            1) redshifts at which output is calculated
            2) with mass as a function of redshift for the main halo
            3) accretion redshifts
            4) subhalo mass at accretion
            5) subhalo mass at z=11 (if 0, the halo collapsed after z=11)
    HISTORY:
       2012-04-16 - Written - Bovy (IAS)
    """
    #First calculate the evolution of the main halo, then calculate the subhalos
    nz= int(math.ceil(zend/dz))+1
    out= []
    out.append(minit)
    m2= minit
    zs= []
    if not _fixdz or _zs is None: z= 0.
    else: z= _zs[0]
    zs.append(z)
    #Set up masses and pre-calculate sig and dsig/dm for speed
    if _logmsgtmres is None:
        logmsgtmres= numpy.linspace(numpy.log(mres),numpy.log(minit),
                                    _NMSINTEGRATE)
        logmsltmres= numpy.linspace(numpy.log(10.**-6.),numpy.log(mres),
                                    _NMSINTEGRATE)
    else:
        logmsgtmres= _logmsgtmres
        logmsltmres= _logmsltmres
    msgtmres= numpy.exp(logmsgtmres)
    msltmres= numpy.exp(logmsltmres)
    if _siggtmres is None:
        if dm is None: #Use standard from cosmolopy
            siggtmres= cp.sigma_r(cp.mass_to_radius(msgtmres,**cosmo),0.,**cosmo)[0]
            sigltmres= cp.sigma_r(cp.mass_to_radius(msltmres,**cosmo),0.,**cosmo)[0]
            dsigdmgtmres= _dsig12dm(msgtmres,siggtmres)
            dsigdmltmres= _dsig12dm(msltmres,sigltmres)
        else:
            #Calculate sigma for gtmres
            thisS= numpy.zeros((len(msgtmres)))
            for ii in range(len(msgtmres)):
                thisS[ii]= massfunction.S(transfer.mToR(msgtmres[ii]),
                                          0.,0.,
                                          window='tophat',
                                          dm=dm,
                                          dmks=dmks)
            #First use thisS for derivative
            dSdm= (thisS-numpy.roll(thisS,-1))/(msgtmres-numpy.roll(msgtmres,-1))
            dSdm[-1]= dSdm[-2] #Hack that does not matter
            dSdm= numpy.fabs(dSdm)
            dsigdmgtmres= dSdm
            siggtmres= numpy.sqrt(thisS)
            #Calculate sigma for ltmres
            thisS= numpy.zeros((len(msgtmres)))
            for ii in range(len(msgtmres)):
                thisS[ii]= massfunction.S(transfer.mToR(msltmres[ii]),
                                          0.,0.,
                                          window='tophat',
                                          dm=dm,
                                          dmks=dmks)
            #First use thisS for derivative
            dSdm= (thisS-numpy.roll(thisS,1))/(msltmres-numpy.roll(msltmres,1))
            dSdm[0]= dSdm[1] #Hack that does not matter
            dSdm= numpy.fabs(dSdm)
            dsigdmltmres= dSdm
            sigltmres= numpy.sqrt(thisS)
    else:
        siggtmres= _siggtmres
        sigltmres= _sigltmres
        dsigdmgtmres= _dsigdmgtmres
        dsigdmltmres= _dsigdmltmres
    dmgtmres= (logmsgtmres[1]-logmsgtmres[0])
    dmltmres= (logmsltmres[1]-logmsltmres[0])
    logmres= math.log(mres)
    log106= math.log(10.**-6.)
    ii= 0
    #Draw collapse redshift for this halo, based on eq. 2.26 of Lacey & Cole 1993
    #print "Calculating collapse redshift ..."
    collapsezs= numpy.linspace(zstart+dz,zend,_NCOLLAPSEZS)
    dc2= 1.675/cp.fgrowth(zstart,cosmo['omega_M_0'])
    dccollapsezs= 1.675/cp.fgrowth(collapsezs,cosmo['omega_M_0'])
    m2indx= int(numpy.floor((numpy.log(m2)-logmres)/dmgtmres))
    if minit/2. >= mres:
        m22indx= int(numpy.floor((numpy.log(m2/2.)-logmres)/dmgtmres))
        pzfgtz= numpy.zeros(_NCOLLAPSEZS)
        for ii in range(_NCOLLAPSEZS):
            pzfgtz[ii]= -integrate.trapz(2.*siggtmres[m22indx:m2indx]*msgtmres[m2indx]/msgtmres[m22indx:m2indx]*(dccollapsezs[ii]-dc2)/(siggtmres[m22indx:m2indx]**2.-siggtmres[m2indx]**2.)**1.5*numpy.exp(-(dccollapsezs[ii]-dc2)**2./2./(siggtmres[m22indx:m2indx]**2.-siggtmres[m2indx]**2.)),
                                         siggtmres[m22indx:m2indx])/numpy.sqrt(2.*numpy.pi)
        #Clean negative probabilities at small z
        kk= 0
        while pzfgtz[kk] < pzfgtz[kk+1] and kk < _NCOLLAPSEZS-2:
            pzfgtz[kk]= -1
            kk+= 1
        pzfgtz[(pzfgtz == -1)]= pzfgtz[kk]
        pzfgtz/= pzfgtz[0]
    else:
        m22indx= int(numpy.floor((numpy.log(m2/2.)-log106)/dmltmres))
        pzfgtz= numpy.zeros(_NCOLLAPSEZS)
        for ii in range(_NCOLLAPSEZS):
            pzfgtz[ii]= -integrate.trapz(2.*siggtmres[0:m2indx]*msgtmres[m2indx]/msgtmres[0:m2indx]*(dccollapsezs[ii]-dc2)/(siggtmres[0:m2indx]**2.-siggtmres[m2indx]**2.)**1.5*numpy.exp(-(dccollapsezs[ii]-dc2)**2./2./(siggtmres[0:m2indx]**2.-siggtmres[m2indx]**2.)),
                                         siggtmres[0:m2indx])/numpy.sqrt(2.*numpy.pi)\
                                         -integrate.trapz(2.*sigltmres[m22indx:_NMSINTEGRATE-1]*msgtmres[m2indx]/msltmres[m22indx:_NMSINTEGRATE-1]*(dccollapsezs[ii]-dc2)/(sigltmres[m22indx:_NMSINTEGRATE-1]**2.-siggtmres[m2indx]**2.)**1.5*numpy.exp(-(dccollapsezs[ii]-dc2)**2./2./(sigltmres[m22indx:_NMSINTEGRATE-1]**2.-siggtmres[m2indx]**2.)),
                                                         sigltmres[m22indx:_NMSINTEGRATE-1])/numpy.sqrt(2.*numpy.pi)
        #Clean negative probabilities at small z
        kk= 0
        while kk < _NCOLLAPSEZS-1 and pzfgtz[kk] < pzfgtz[kk+1]:
            pzfgtz[kk]= -1
            kk+= 1
        pzfgtz[(pzfgtz == -1)]= pzfgtz[kk]
        pzfgtz/= pzfgtz[0]
    #return pzfgtz
    #Draw collapse redshift
    kk= 0
    r= numpy.random.uniform()
    while r < pzfgtz[kk] and kk < _NCOLLAPSEZS-1: kk+= 1
    zc= collapsezs[kk]
    #print numpy.log10(m2), zc
    if justcollapsez: return (None,[0.],None,None,None,zc,None)
    #Set up subhalos
    zacc= []
    newhalos= []
    if not _nosubstructure: print "Working on main halo ..."
    foundhalf= False
    while z < zend:
        if _fixdz and not _zs is None: dz= _zs[ii+1]-_zs[ii]
        if z < zstart:
            out.append(-1)
            ii+= 1
            z+= dz
            zs.append(z)
            continue
        #elif z > zc:
        #    out.append(0.)
        #    ii+= 1
        #    z+= dz
        #    zs.append(z)
        #    continue
        #Find closest index for m2/2.
        m2indx= int(numpy.floor((numpy.log(m2)-logmres)/dmgtmres))
        m22indx= int(numpy.floor((numpy.log(m2/2.)-logmres)/dmgtmres))
        if not _nosubstructure:
            #Calculate dndm
            dndlogms= _dndlogm_fast(msgtmres[0:m22indx+1],
                                    msgtmres[m2indx],
                                    siggtmres[0:m22indx+1],
                                    siggtmres[m2indx],
                                    dsigdmgtmres[0:m22indx+1],
                                    z,dz)
            p= integrate.trapz(dndlogms,logmsgtmres[0:m22indx+1],dx=dmgtmres)
            r= numpy.random.uniform()
            if r > p: #Do not fragment
                newmass= 0.
            else: #Fragment
                zerolist= [0.]
                zerolist.extend(list(integrate.cumtrapz(dndlogms,
                                                        logmsgtmres[0:m22indx+1],
                                                        dx=dmgtmres)))
                cumulmass= numpy.array(zerolist)
                cumulmass/= p #normalize
            #Now draw a mass
                kk= 0
                r= numpy.random.uniform()
                while r > cumulmass[kk] and kk < m22indx: kk+= 1
                newmass= msgtmres[kk]
                zacc.append(z)
                newhalos.append(newmass)
        #Calculate accretion fraction
        if not _nosubstructure:
            dndlogms= _dndlogm_fast(msltmres,msgtmres[m2indx],
                                    sigltmres,siggtmres[m2indx],
                                    dsigdmltmres,z,dz)
            f= integrate.trapz(dndlogms*msltmres,logmsltmres,dx=dmltmres)/m2
        else:
            if m22indx < 0.:
                m22indx= int(numpy.floor((numpy.log(m2/2.)-log106)/dmltmres))
                m22ltzero= True
            else:
                m22ltzero= False
            if m2indx < 0.:
                m2indx= int(numpy.floor((numpy.log(m2)-log106)/dmltmres))
                sigm2= sigltmres[m2indx]
            else:
                sigm2= siggtmres[m2indx]
            if m22ltzero:
                dndlogms= _dndlogm_fast(msltmres[0:m22indx+1],
                                        m2,
                                        sigltmres[0:m22indx+1],
                                        sigm2,
                                        dsigdmltmres[0:m22indx+1],
                                        z,dz)
                f= integrate.trapz(dndlogms*msltmres[0:m22indx+1],
                                   logmsltmres[0:m22indx+1],
                                   dx=dmltmres)/m2
            else:
                #Calculate dndm
                dndlogms= _dndlogm_fast(msltmres,msgtmres[m2indx],
                                        sigltmres,siggtmres[m2indx],
                                        dsigdmltmres,z,dz)
                f= integrate.trapz(dndlogms*msltmres,logmsltmres,dx=dmltmres)/m2
                #Calculate dndm
                dndlogms= _dndlogm_fast(msgtmres[0:m22indx+1],
                                        msgtmres[m2indx],
                                        siggtmres[0:m22indx+1],
                                        siggtmres[m2indx],
                                        dsigdmgtmres[0:m22indx+1],
                                        z,dz)
                f+= integrate.trapz(dndlogms*msgtmres[0:m22indx+1],
                                    logmsgtmres[0:m22indx+1],
                                    dx=dmgtmres)/m2
        if _nosubstructure: newmass= 0. #HACK
        m2= m2*(1.-f)-newmass
        #if ii%100 == 0 and not newmass == 0.: print z, _nosubstructure*p, numpy.log10(newmass), numpy.log10(f), numpy.log10(m2)
        #elif ii%100 == 0 and not _nosubstructure: print z, p, numpy.log10(f), numpy.log10(m2)
        #elif ii%100 == 0: print z, numpy.log10(f), numpy.log10(m2)
        #if ii%100 == 0:
        #    sys.stdout.write('\r'+"z = %.4f ..." % z)
        #    sys.stdout.flush()
        #Formation redshift?
        #if m2 < minit/2. and not foundhalf: 
        #    zc= z-dz
        #    print numpy.log10(minit), zc
        #    foundhalf= True
        out.append(m2)
        ii+= 1
        z+= dz
        zs.append(z)
        #Larger or smaller dz?
        if not _nosubstructure and p < 0.1 and not _fixdz: dz*= 2.
        if not _nosubstructure and p > 0.3 and not _fixdz: dz/= 2.
    #sys.stdout.write('\r'+_ERASESTR+'\r')
    #sys.stdout.flush()
    newhalos= numpy.array(newhalos)
    #Now work on the subhalos
    if not _nosubstructure:
        nnew= len(zacc)
        newhalosprerei= numpy.zeros(nnew)
        subzcs= numpy.zeros(nnew)
        for ii in range(nnew):
            if dontevolvesub:
                newhalosprerei[ii]= newhalos[ii]
                subzcs[ii]= -1.
                continue
            if not _nosubstructure:
                sys.stdout.write('\r'+"Working on subhalo %i / %i ..." % (ii+1,nnew))
                sys.stdout.flush()
            #Recursion
            dummyzs,out, dummyzacc,dummynewhalos, dummynewhalosprerei,subzc,dummyzc= mergertree(minit=newhalos[ii],
                                                                                                zend=zend,
                                                                                                mres=mres,
                                                                                                dz=dzsub,
                                                                                                zstart=zacc[ii],
                                                                                                justcollapsez=False,
                                                                                                _nosubstructure=True,
                                                                                                _fixdz=True,
                                                                                                _logmsgtmres=logmsgtmres,
                                                                                                _logmsltmres=logmsltmres,
                                                                                                
                                                                                                _siggtmres=siggtmres,
                                                                                                _sigltmres=sigltmres,
                                                                                                _dsigdmgtmres=dsigdmgtmres,
                                                                                                _dsigdmltmres=dsigdmltmres)
            newhalosprerei[ii]= out[-1]
            subzcs[ii]= subzc
        if not _nosubstructure:
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
        #Deal with ahlos accreted right before end
        indx= (newhalosprerei == -1)
        newhalosprerei[indx]= newhalos[indx]
    else:
        newhalosprerei= None
        subzcs= None
    return (zs,out,zacc,newhalos,newhalosprerei,zc,subzcs)

def _epsP(m2,mres,z,dz):
    ms= numpy.linspace(numpy.log(mres),numpy.log(m2/2.),_NMSINTEGRATE)
    dndlogms= _dndlogm(numpy.exp(ms),m2,z,dz)
    return integrate.trapz(dndlogms,ms,dx=(ms[1]-ms[0]))

def _epsF(m2,mres,z,dz):
    ms= numpy.linspace(numpy.log(10.**-6.),numpy.log(mres),_NMSINTEGRATE)
    dndlogms= _dndlogm(numpy.exp(ms),m2,z,dz)*numpy.exp(ms)/m2
    return integrate.trapz(dndlogms,ms,dx=(ms[1]-ms[0]))

def _dndlogm_fast(m1,m2,sig1,sig2,dsig1dm,z,dz):
    return numpy.fabs(1./math.sqrt(2.*math.pi)/(sig1**2.-sig2**2.)**1.5\
                          *_ddcdz(z,dz)\
                          *dsig1dm\
                          *dz\
                          *m2)

def _dndlogm(m1,m2,z,dz):
    sig1= cp.sigma_r(cp.mass_to_radius(m1,**cosmo),0.,**cosmo)[0]
    sig2= cp.sigma_r(cp.mass_to_radius(m2,**cosmo),0.,**cosmo)[0]
    return numpy.fabs(1./math.sqrt(2.*math.pi)/(sig1**2.-sig2**2.)**1.5\
                          *_ddcdz(z,dz)\
                          *_dsig12dm(m1,sig1)\
                          *dz\
                          *m2)

def _ddcdz(z,dz):
    dc1= 1./cp.fgrowth(z+dz,cosmo['omega_M_0'])
    dc2= 1./cp.fgrowth(z,cosmo['omega_M_0'])
    return 1.675*(dc1-dc2)/dz

def _dsig12dm(m,sig):
    sig2= cp.sigma_r(cp.mass_to_radius(m*(1+10.**-6.),**cosmo),0.,**cosmo)[0]
    return 2.*sig*(sig-sig2)/m/10.**-6. #dm=1 solar mass
