import sys
import os, os.path
import copy
from optparse import OptionParser
import cPickle as pickle
import numpy
from scipy import interpolate, integrate
try:
    from galpy.util import bovy_plot
except ImportError:
    import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import transfer
_NEW= True
def nfw(r,rs=30.):
    return 1./r/(1.+r/rs)**2.
def completenessIntegral(r,phi,theta,rvir,fsdss,dmax2):
    #Calculate whether d < dmax
    x= r*numpy.cos(theta)*numpy.sin(phi)
    y= r*numpy.sin(theta)*numpy.sin(phi)
    z= r*numpy.cos(phi)
    d2= (x-8.)**2.+y**2.+z**2.
    #Smooth for integration purposes
    return 1./(1.+numpy.exp((d2-dmax2)/1.))*r**2.*numpy.sin(phi)*nfw(r)*r #Anti-bias, Diemand et al. (2007)
    #if d2 < dmax2: return r**2.*numpy.sin(phi)*nfw(r)*r
    #else: return 0.

def sdss_completeness(mv,rvir=261.,fsdss=0.194,koposov=True):
    """SDSSS completeness, rvir in kpc"""
    #Calculate rmax
    if koposov:
        dmax= 10.**(1.1-0.228*mv) #kpc
    else:
        dmax= (3./4./numpy.pi/fsdss)**(1./3.)*10.**((-0.6*mv-5.23)/3.+3.) #kpc
    """
    #Rough approximation
    if isinstance(mv,float):
        if dmax > rvir: return fsdss
        else: return fsdss*(dmax/rvir)**3.
    elif isinstance(mv,numpy.ndarray):
        out= numpy.zeros_like(mv)
        out[:]= fsdss
        indx= (dmax <= rvir)
        out[indx]= fsdss*(dmax[indx]/rvir)**3.
        return out
    """
    if dmax > 100.: #Approximation that seems okay
        return fsdss*integrate.tplquad(completenessIntegral,0.,2.*numpy.pi,
                                       lambda x: 0., lambda x: numpy.pi,
                                       lambda x,y: 0., 
                                       lambda x,y: numpy.amin([dmax,rvir]),
                                       (rvir,fsdss,1000.**2.))[0]/\
                                       integrate.tplquad(completenessIntegral,
                                                         0.,2.*numpy.pi,
                                                         lambda x: 0., 
                                                         lambda x: numpy.pi,
                                                         lambda x,y: 0., 
                                                         lambda x,y: rvir,
                                                         (rvir,fsdss,1000.**2.))[0]
    else:
        return fsdss*integrate.tplquad(completenessIntegral,0.,2.*numpy.pi,
                                       lambda x: 0., lambda x: numpy.pi,
                                       lambda x,y: 0., lambda x,y: rvir,
                                       (rvir,fsdss,dmax**2.))[0]/\
                                       integrate.tplquad(completenessIntegral,
                                                         0.,2.*numpy.pi,
                                                         lambda x: 0., 
                                                         lambda x: numpy.pi,
                                                         lambda x,y: 0., 
                                                         lambda x,y: rvir,
                                                         (rvir,fsdss,1000.**2.))[0]
                                                        
def plot_mergertree(parser):
    (options,args)= parser.parse_args()
    if len(args) < 1 or options.plotfile is None:
        parser.print_help()
        return
    #Restore all files
    nargs= len(args)
    outs= []
    zaccs= []
    newhaloss= []
    newhalosprereis= []
    zcmains= []
    subzcs= []
    for ii in range(nargs):
        savefile= open(args[ii],'rb')
        zs= pickle.load(savefile)
        out= pickle.load(savefile)
        zacc= pickle.load(savefile)
        newhalos= pickle.load(savefile)
        newhalosprerei= pickle.load(savefile)
        zcmain= pickle.load(savefile)
        subzc= pickle.load(savefile)
        savefile.close()
        outs.append(out)
        zaccs.append(zacc)
        newhaloss.append(newhalos)
        newhalosprereis.append(newhalosprerei)
        zcmains.append(zcmain)
        subzcs.append(subzc)
    #Now plot
    if options.plottype.lower() == 'zacc':
        nbins= 31
        bovy_plot.bovy_print()
        bovy_plot.bovy_hist(zaccs[0],range=[0.,11.],
                            xlabel=r'$z_{\mathrm{acc}}$',
                            bins=nbins,color='k',histtype='step')
        if nargs > 1:
            for ii in range(1,nargs):
                bovy_plot.bovy_hist(zaccs[ii],range=[0.,11.],bins=nbins,
                                    overplot=True,
                                    color='0.65',histtype='step')
        _ADDCUMUL= True
        if _ADDCUMUL:
            ax= pyplot.gca()
            ylimits= ax.get_ylim()
            cumulzs= numpy.ones_like(zaccs[0])
            bovy_plot.bovy_plot(zaccs[0],
                                numpy.cumsum(cumulzs)/len(zaccs[0])*ylimits[1],
                                overplot=True,
                                color='k',ls='-')
            if nargs > 1:
                for ii in range(1,nargs):
                    cumulzs= numpy.ones_like(zaccs[ii])
                    bovy_plot.bovy_plot(zaccs[ii],
                                        numpy.cumsum(cumulzs)/len(zaccs[ii])*ylimits[1],
                                        overplot=True,
                                        color='0.65',ls='-')
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)
    elif options.plottype.lower() == 'mass':
        nbins= 61
        #Calculate mass function
        ii= 0
        lnms= numpy.log(newhaloss[ii])
        hist, edges= numpy.histogram(lnms,
                                     range=[numpy.log(10.**5.),
                                            numpy.log(10.**11.)],
                                     bins=nbins)
        hist= hist.astype('float')
        xs= (edges+numpy.roll(edges,-1))/2.
        hist/= numpy.fabs(xs[1]-xs[0])
        if not options.saveavg is None:
            allhists= []
            allhists.append(copy.copy(hist))
        if options.addavg:
            thiscolor= '0.65'
            avghist= hist
        else:
            thiscolor= 'k'
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(numpy.exp(xs[0:-1]),hist,'-',loglog=True,
                            color=thiscolor,
                            xrange=[10.**6.,10.**11.],
                            yrange=[0.1,10.**4.],
                            xlabel=r'$M_{\mathrm{halo}} / M_\odot$',
                            ylabel=r'$\mathrm{d} N / \mathrm{d} \ln M_{\mathrm{halo}}$')
        if nargs > 1:
            for ii in range(1,nargs):
                lnms= numpy.log(newhaloss[ii])
                hist, edges= numpy.histogram(lnms,
                                             range=[numpy.log(10.**5.),
                                                    numpy.log(10.**11.)],
                                             bins=nbins)
                hist= hist.astype('float')
                xs= (edges+numpy.roll(edges,-1))/2.
                hist/= numpy.fabs(xs[1]-xs[0])
                if options.addavg:
                    avghist+= hist
                if not options.saveavg is None:
                    allhists.append(copy.copy(hist))
                bovy_plot.bovy_print()
                bovy_plot.bovy_plot(numpy.exp(xs[0:-1]),hist,'-',loglog=True,
                                    overplot=True,color='0.65')
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        if options.addavg:
            bovy_plot.bovy_plot(numpy.exp(xs[0:-1]),
                                avghist/nargs,'-',loglog=True,
                                overplot=True,color='k',lw=1.5)
        bovy_plot.bovy_end_print(options.plotfile)
        if options.addavg and not options.saveavg is None:
            #Calculate std
            dev= numpy.zeros((2,len(avghist)))
            for ii in range(len(avghist)):
                thish= [h[ii] for h in allhists]
                #print len(thish), avghist[ii]/nargs, thish
                #dev[ii]= numpy.std(thish)
                dev[0,ii]= sorted(thish)[int(numpy.round(0.16*nargs))]
                dev[1,ii]= sorted(thish)[int(numpy.round(0.84*nargs))]
            savefile= open(options.saveavg,'wb')
            pickle.dump(xs[0:-1],savefile)
            pickle.dump(avghist/nargs,savefile)
            pickle.dump(dev,savefile)
            savefile.close()
    elif options.plottype.lower() == 'z11mass':
        nbins= 61
        #Calculate mass function
        ii= 0
        lnms= numpy.log(newhalosprereis[ii])
        hist, edges= numpy.histogram(lnms,
                                     range=[numpy.log(10.**5.),
                                            numpy.log(10.**11.)],
                                     bins=nbins) 
        hist= hist.astype('float')
        xs= (edges+numpy.roll(edges,-1))/2.
        hist/= numpy.fabs(xs[1]-xs[0])
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(numpy.exp(xs[0:-1]),hist,'k-',loglog=True,
                            xrange=[10.**5.,10.**11.],
                            yrange=[0.1,10.**4.],
                            xlabel=r'$M_{\mathrm{halo}} / M_\odot$',
                            ylabel=r'$\mathrm{d} N / \mathrm{d} \ln M_{\mathrm{halo}}$')
        if nargs > 1:
            for ii in range(1,nargs):
                lnms= numpy.log(newhalosprereis[ii])
                hist, edges= numpy.histogram(lnms,
                                             range=[numpy.log(10.**5.),
                                                    numpy.log(10.**11.)],
                                             bins=nbins)
                hist= hist.astype('float')
                xs= (edges+numpy.roll(edges,-1))/2.
                hist/= numpy.fabs(xs[1]-xs[0])
                bovy_plot.bovy_print()
                bovy_plot.bovy_plot(numpy.exp(xs[0:-1]),hist,'-',loglog=True,
                                    overplot=True,color='0.65')
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)
    elif options.plottype.lower() == 'z0massvsz11mass':
        nbins= 61
        #For each tree, calculate the mean z11mass
        z0masses= numpy.exp(numpy.linspace(numpy.log(10.**6.),numpy.log(10.**11.),nbins+1))
        thismean= numpy.zeros(nbins)
        ii= 0
        for jj in range(nbins):
            thisindx= (newhaloss[ii] > z0masses[jj])*(newhaloss[ii] <= z0masses[jj+1])
            thismean[jj]= numpy.mean(newhalosprereis[ii][thisindx])
        if options.addavg:
            thiscolor= '0.65'
            avgmean= thismean
        else:
            thiscolor= 'k'
        xs= (z0masses+numpy.roll(z0masses,-1))/2.
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(xs[0:-1],thismean,'-',loglog=True,
                            color=thiscolor,
                            xrange=[10.**4.5,10.**11.],
                            yrange=[10.**4.5,10.**11.],
                            xlabel=r'$M_{\mathrm{halo}}(z=z_{\mathrm{acc}}) / M_\odot$',
                            ylabel=r'$M_{\mathrm{halo}}(z=11) / M_\odot$')
        if nargs > 1:
            for ii in range(1,nargs):
                for jj in range(nbins):
                    thisindx= (newhaloss[ii] > z0masses[jj])*(newhaloss[ii] <= z0masses[jj+1])
                    thismean[jj]= numpy.mean(newhalosprereis[ii][thisindx])
                if options.addavg:
                    avgmean+= thismean
                bovy_plot.bovy_plot(xs[0:-1],thismean,'-',color='0.65',
                                    overplot=True,loglog=True)
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        #Add average
        if options.addavg:
            bovy_plot.bovy_plot(xs[0:-1],avgmean/nargs,'k-',lw=1.5,
                                overplot=True,loglog=True)
        #Overlay a few lines
        bovy_plot.bovy_plot([10.**4.5,10.**11.],[10.**4.5,10.**11.],'k-',lw=0.5,
                            overplot=True)
        bovy_plot.bovy_plot([10.**4.5,10.**11.],[10.**3.5,10.**10.],'k-',lw=0.5,
                            overplot=True)
        bovy_plot.bovy_plot([10.**4.5,10.**11.],[10.**5.5,10.**12.],'k-',lw=0.5,
                            overplot=True)
        bovy_plot.bovy_plot([10.**4.5,10.**11.],[10.**2.5,10.**9.],'k-',lw=0.5,
                            overplot=True)
        bovy_plot.bovy_plot([10.**4.5,10.**11.],[10.**6.5,10.**13.],'k-',lw=0.5,
                            overplot=True)
        bovy_plot.bovy_end_print(options.plotfile)
    elif (options.plottype.lower() == 'lumfunc' or options.plottype.lower() == 'sdsslumfunc') and not options.cumul:
        nbins= 61
        mvmassnbins= 31
        #Calculate ms for each halo
        ii= 0
        mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                 zaccs[ii],
                                                 options)
        #Apply fs
        #mstarspostrei*= options.fs/0.01
        mstarsprerei*= options.fs/0.01
        #Calculate MV
        mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                               mvpersolarmass=4.8)
        mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
#                              mvpersolarmass=4.8)
        #Combine the two
        mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
        #Calculate histogram
        hist, edges= numpy.histogram(mv,
                                     range=[-10.,0.],
                                     bins=nbins)
        hist= hist.astype('float')
        xs= (edges+numpy.roll(edges,-1))/2.
        #Calculate the mv-mhalo relationship
        mvs= numpy.linspace(-10.,0.,mvmassnbins+1)
        mvmass= [[] for jj in range(mvmassnbins)]
        for jj in range(mvmassnbins):
            thisindx= (mv > mvs[jj])*(mv <= mvs[jj+1])
            mvmass[jj].extend(newhaloss[ii][thisindx])
        if options.plottype.lower() == 'sdsslumfunc':
            #Get or calculate the selection fraction
            if not options.sdssselect is None and os.path.exists(options.sdssselect):
                savefile= open(options.sdssselect,'rb')
                sdssselect= pickle.load(savefile)
                savefile.close()
            else:
                #Calculate!
                sdssselect= numpy.zeros(len(xs)-1)
                for ii in range(len(xs)-1):
                    sdssselect[ii]= sdss_completeness(xs[ii])
                if not options.sdssselect is None:
                    savefile= open(options.sdssselect,'wb')
                    pickle.dump(sdssselect,savefile)
                    pickle.dump(xs,savefile)
                    savefile.close()
            hist*= sdssselect
            yrange= [0.1,100.]
        else:
            yrange= [1.,10000.]
        hist/= numpy.fabs(xs[1]-xs[0])
        if not options.saveavg is None:
            allhists= []
            allhists.append(copy.copy(hist))
        bovy_plot.bovy_print()
        if options.addavg:
            thiscolor= '0.65'
            avghist= hist
        else:
            thiscolor= 'k'
        bovy_plot.bovy_plot(xs[0:-1],hist,'-',semilogy=True,
                            color=thiscolor,
                            xrange=[0.,-10.],
                            yrange=yrange,
                            xlabel=r'$M_V$',
                            ylabel=r'$\mathrm{d} N / \mathrm{d} M_V$')
        if nargs > 1:
            for ii in range(1,nargs):
                mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                         zaccs[ii],
                                                         options)
                #Apply fs
                #mstarspostrei*= options.fs/0.01
                mstarsprerei*= options.fs/0.01
                #Calculate MV
                mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                                       mvpersolarmass=4.8)
                mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
#                                      mvpersolarmass=4.8)
                #Combine the two
                mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
                #Calculate histogram
                hist, edges= numpy.histogram(mv,
                                             range=[-10.,0.],
                                             bins=nbins)
                hist= hist.astype('float')
                xs= (edges+numpy.roll(edges,-1))/2.
                #Calculate Mv-mass relation
                for jj in range(mvmassnbins):
                    thisindx= (mv > mvs[jj])*(mv <= mvs[jj+1])
                    mvmass[jj].extend(newhaloss[ii][thisindx])
                if options.plottype.lower() == 'sdsslumfunc':
                    hist*= sdssselect
                hist/= numpy.fabs(xs[1]-xs[0])
                if options.addavg:
                    avghist+= hist
                if not options.saveavg is None:
                    allhists.append(copy.copy(hist))
                bovy_plot.bovy_plot(xs[0:-1],hist,'-',semilogy=True,
                                    overplot=True,color='0.65')
        #Add average
        if options.addavg:
            bovy_plot.bovy_plot(xs[0:-1],avghist/nargs,'-',semilogy=True,
                                overplot=True,color='k',lw=1.5)
            if not options.saveavg is None:
                #Calculate std
                dev= numpy.zeros((2,len(avghist)))
                for ii in range(len(avghist)):
                    thish= [h[ii] for h in allhists]
                    #print len(thish), avghist[ii]/nargs, thish
                    #dev[ii]= numpy.std(thish)
                    dev[0,ii]= sorted(thish)[int(numpy.round(0.16*nargs))]
                    dev[1,ii]= sorted(thish)[int(numpy.round(0.84*nargs))]
                savefile= open(options.saveavg,'wb')
                pickle.dump(xs[0:-1],savefile)
                pickle.dump(avghist/nargs,savefile)
                pickle.dump(dev,savefile)
                pickle.dump(allhists,savefile)
                pickle.dump(mvs,savefile)
                pickle.dump([numpy.median(mvmass[jj]) for jj in range(mvmassnbins)],savefile)
                savefile.close()
        #Add Koposov points
        if options.plottype.lower() == 'sdsslumfunc':
            koposov= numpy.loadtxt('koposov09a.txt',delimiter='|')
            bovy_plot.bovy_plot(koposov[:,0],koposov[:,1],'Dw',overplot=True,
                                mew=1.)
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)
    elif options.plottype.lower() == 'lumfunc' or options.plottype.lower() == 'sdsslumfunc' and options.cumul:
        #Calculate ms for each halo
        ii= 0
        mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                 zaccs[ii],
                                                 options)
        #Apply fs
        #mstarspostrei*= options.fs/0.01
        mstarsprerei*= options.fs/0.01
        #Calculate MV
        mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                               mvpersolarmass=4.8)
        mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
        #Combine the two
        mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
        #Calculate cumulative
        mv= sorted(mv)
        #Load sdssselect
        if options.plottype.lower() == 'sdsslumfunc':
            #Get the selection fraction
            if os.path.exists(options.sdssselect):
                savefile= open(options.sdssselect,'rb')
                sdssselect= pickle.load(savefile)
                sdssselectxs= pickle.load(savefile)
                savefile.close()
                selectInterp= interpolate.InterpolatedUnivariateSpline(sdssselectxs[0:-1],numpy.log(sdssselect))
            else:
                raise IOError("sdssselect file does not exist, calculate first with plottype=lumfunc")
        #Calculate sdssselect for each mv
        if options.plottype.lower() == 'sdsslumfunc':
            select= numpy.exp(selectInterp(mv))
            yrange= [1.,300.]
        else:
            select= numpy.ones_like(mv)
            yrange= [1.,1000.]
        cumuldist= numpy.cumsum(select)
        bovy_plot.bovy_print()
        if options.addavg:
            thiscolor= '0.65'
            allmv= []
            allselect= []
            allmv.extend(mv)
            allselect.extend(select)
        else:
            thiscolor= 'k'
        bovy_plot.bovy_plot(mv,cumuldist,'-',semilogy=True,
                            color=thiscolor,
                            xrange=[0.,-10.],
                            yrange=yrange,
                            xlabel=r'$M_V$',
                            ylabel=r'$N(<M_V)$')
        if nargs > 1:
            for ii in range(1,nargs):
                mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                         zaccs[ii],
                                                         options)
                #Apply fs
                #mstarspostrei*= options.fs/0.01
                mstarsprerei*= options.fs/0.01
                #Calculate MV
                mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                                       mvpersolarmass=4.8)
                mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
                #Combine the two
                mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
                mv= sorted(mv)
                #Calculate sdssselect for each mv
                if options.plottype.lower() == 'sdsslumfunc':
                    select= numpy.exp(selectInterp(mv))
                else:
                    select= numpy.ones_like(mv)
                cumuldist= numpy.cumsum(select)
                if options.addavg:
                    allmv.extend(mv)
                    allselect.extend(select)
                bovy_plot.bovy_plot(mv,cumuldist,'-',semilogy=True,
                                    color='0.65',overplot=True)
        #Add average
        if options.addavg:
            sortindx= numpy.argsort(numpy.array(allmv))
            mv= sorted(allmv)
            select= numpy.array(allselect).flatten()[sortindx]
            cumuldist= numpy.cumsum(select)
            bovy_plot.bovy_plot(mv,cumuldist/nargs,'-',semilogy=True,
                                overplot=True,color='k',lw=1.5)
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)
    elif options.plottype.lower() == 'mvmass':
        nbins= 61
        #Calculate ms for each halo
        ii= 0
        mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                 zaccs[ii],
                                                 options)
        #Apply fs
        #mstarspostrei*= options.fs/0.01
        mstarsprerei*= options.fs/0.01
        #Calculate MV
        mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                               mvpersolarmass=4.8)
        mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
#                              mvpersolarmass=4.8)
        #Combine the two
        mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
        #Calculate the mv-mhalo relationship
        mvs= numpy.linspace(-10.,0.,nbins+1)
        thismean= numpy.zeros(nbins)
        for jj in range(nbins):
            thisindx= (mv > mvs[jj])*(mv <= mvs[jj+1])
            thismean[jj]= numpy.mean(newhaloss[ii][thisindx])
        if options.addavg:
            thiscolor= '0.65'
            avgmean= thismean
        else:
            thiscolor= 'k'
        xs= (mvs+numpy.roll(mvs,-1))/2.
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(xs[0:-1],thismean,'-',semilogy=True,
                            color=thiscolor,
                            xrange=[0.,-10.],
                            yrange=[10.**6.,10.**11.],
                            xlabel=r'$M_V$',
                            ylabel=r'$M_{\mathrm{halo}} / M_\odot$')
        if nargs > 1:
            for ii in range(1,nargs):
                mstarspostrei,mstarsprerei= calc_z0mstar(newhaloss[ii],newhalosprereis[ii],
                                                         zaccs[ii],
                                                         options)
                #Apply fs
                #mstarspostrei*= options.fs/0.01
                mstarsprerei*= options.fs/0.01
                #Calculate MV
                mvpostrei= transfer.mv(None,None,None,thismstar=mstarspostrei,
                                       mvpersolarmass=4.8)
                mvprerei= transfer.mv(None,None,None,thismstar=mstarsprerei)
#                                      mvpersolarmass=4.8)
                #Combine the two
                mv= -2.5*numpy.log10(10.**(-mvpostrei/2.5)+10.**(-mvprerei/2.5))
                #Calculate relation
                for jj in range(nbins):
                    thisindx= (mv > mvs[jj])*(mv <= mvs[jj+1])
                    thismean[jj]= numpy.mean(newhaloss[ii][thisindx])
                if options.addavg:
                    avgmean+= thismean
                bovy_plot.bovy_plot(xs[0:-1],thismean,'-',semilogy=True,
                                    overplot=True,color='0.65')
        #Add vbc label
        if options.vbc == 0:
            bovy_plot.bovy_text(r'$v_{bc} = 0$',top_right=True,size=14.)
        else:
            bovy_plot.bovy_text(r'$v_{bc} = %i\,\sigma_{bc}$' % options.vbc,
                                top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)
    return None

def calc_z0mstar(mass,massprerei,zacc,options):
    """Returns stellar mass (after,before) reionization"""
    masspostrei= mass-massprerei
    #Load pre-rei Mstar/Mhalo relation
    if os.path.exists(options.masssav):
        #Saved, restore
        savefile= open(options.masssav,'rb')
        ms= pickle.load(savefile)
        dmks= pickle.load(savefile)
        dndmv0= pickle.load(savefile)
        dm0= pickle.load(savefile)
        dndmv1= pickle.load(savefile)
        dm1= pickle.load(savefile)
        dndmv2= pickle.load(savefile)
        dm2= pickle.load(savefile)
        dndmv4= pickle.load(savefile)
        dm4= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("masssav does not exist; calculate using plot_massfunction.py first")
    if os.path.exists(options.lumfuncsav):
        savefile= open(options.lumfuncsav,'rb')
        mstar0= pickle.load(savefile)
        mv0= pickle.load(savefile)
        dmstardm0= pickle.load(savefile)
        dmvdmstar0= pickle.load(savefile)
        mf0= pickle.load(savefile)
        rlss0= pickle.load(savefile)
        mstar1= pickle.load(savefile)
        mv1= pickle.load(savefile)
        dmstardm1= pickle.load(savefile)
        dmvdmstar1= pickle.load(savefile)
        mf1= pickle.load(savefile)
        rlss1= pickle.load(savefile)
        mstar2= pickle.load(savefile)
        mv2= pickle.load(savefile)
        dmstardm2= pickle.load(savefile)
        dmvdmstar2= pickle.load(savefile)
        mf2= pickle.load(savefile)
        rlss2= pickle.load(savefile)
        mstar4= pickle.load(savefile)
        mv4= pickle.load(savefile)
        dmstardm4= pickle.load(savefile)
        dmvdmstar4= pickle.load(savefile)
        mf4= pickle.load(savefile)
        rlss4= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("lumfuncsav not given, need to calculate this with plot_lumfunc first")
    #Load correct vbc
    if options.vbc == 0:
        mstar= mstar0
        mf= mf0
        rlss= rlss0
    elif options.vbc == 1:
        mstar= mstar1
        mf= mf1
        rlss= rlss1
    elif options.vbc == 2:
        mstar= mstar2
        mf= mf2
        rlss= rlss2
    elif options.vbc == 3:
        mstar= mstar4
        mf= mf4
        rlss= rlss4
    #Interpolate this relation
    #indx= (ms < 10.**13.)
    #ms= ms[indx]
    #mstar= mstar[indx]
    #highmassratio= mstar[numpy.argmax(ms)]/ms[numpy.argmax(ms)]
    #mstarmsInterp= interpolate.interp1d(numpy.log(ms)[::-1],
    #                                    numpy.log(mstar/ms/highmassratio)[::-1],
    #                                    kind='linear',bounds_error=False)
    #out= massprerei*numpy.exp(mstarmsInterp(numpy.log(massprerei)))*highmassratio
    out= transfer.mstar(massprerei,None,None,
                        thismf=mf,rlss=rlss,new=_NEW)
    #Add z < zrei contribution
    vcircs= transfer.vcirc(mass,numpy.array(zacc))
    return (0.001/6.25*masspostrei/(1./+0.26*(options.vcrit/vcircs)**3.)**3.,
            out)
    
def get_options():
    usage = "usage: %prog [options] <savefilenames>\n\nsavefilenames= name of the files that stuff was be saved to in run_mergertree.py (as many as you want)"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    parser.add_option("-t","--plottype",dest='plottype',default='mass',
                      help="Type of plot to make ('mass' for z=0 mass function, 'z11mass' for z=11 mass function, 'lumfunc' for z=0 luminosity function, 'zacc' for histogram of accretion redshifts)")
    parser.add_option("--masssav",dest='masssav',default=None,
                      help="name of the file that holds the mass function pickles")
    parser.add_option("--lumfuncsav",dest='lumfuncsav',default=None,
                      help="name of the file that holds the luminosity function pickles")
    parser.add_option("--vbc",dest='vbc',default=0,type='int',
                      help="Vbc in sigma (0,1,2,3)")
    parser.add_option("--fs",dest='fs',default=.03,type='float',
                      help="Star formation efficiency")
    parser.add_option("--vcrit",dest='vcrit',default=35.,type='float',
                      help="Critical circular velocity for star formation after reionization")
    parser.add_option("--addavg",action="store_true", dest="addavg",
                      default=False,
                      help="Add the average effect")
    parser.add_option("--cumul",action="store_true", dest="cumul",
                      default=False,
                      help="Plot the cumulative distribution, for lumfunc or sdsslumfunc")
    parser.add_option("--sdssselect",dest='sdssselect',default=None,
                      help="Save the SDSS seletion function to this file")
    parser.add_option("--saveavg",dest='saveavg',default=None,
                      help="If set, save the average to this file")
    return parser

if __name__ == '__main__':
    plot_mergertree(get_options())
    
