import sys
import os, os.path
from optparse import OptionParser
import cPickle as pickle
import numpy
from scipy import interpolate
import transfer
import massfunction
try:
    from galpy.util import bovy_plot
except ImportError:
    import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from matplotlib.patches import FancyArrowPatch
from plot_massfunction import calc_dndm
_ERASESTR= "                                                                                "
_NEW= True
def calc_mstaretc(ms,options,vbc):
    if vbc == 0.: ncost= 1
    else: ncost= options.ncost
    #First calculate thismf and rlss
    thismf, rlss= transfer.mf(options.zrei,vbc,ncost=ncost,returnrlss=True,
                              new=_NEW)
    mstar= transfer.mstar(ms,options.zrei,vbc,ncost=ncost,
                          thismf=thismf,rlss=rlss,new=_NEW)
    mv= transfer.mv(ms,options.zrei,
                    vbc,ncost=ncost,thismstar=mstar)
    dmstardm= massfunction.dmstardm(ms,options.zrei,vbc,ncost=ncost,
                                    thismf=thismf,rlss=rlss)
    dmvdmstar= massfunction.dmvdmstar(mstar)
    return mstar, mv, dmstardm, dmvdmstar, thismf, rlss

def pvbc(vbc):
    return (3./2./numpy.pi/transfer.sbc**2.)**1.5*4.*numpy.pi*vbc**2.*numpy.exp(-1.5*vbc**2./transfer.sbc**2.)

def calc_avg(ms,options):
    vbcs= numpy.linspace(0.,2.5,options.nvbcavg)*transfer.sbc
    mvs= numpy.linspace(5.,-15.,1001)
    avg= numpy.zeros_like(mvs)
    norm= 0.
    for ii in range(1,options.nvbcavg):
        print "Working on %i / %i ..." % (ii+1,options.nvbcavg)
        mstar, mv, dmstardm, dmvdmstar,dummymf,dummyrlss= calc_mstaretc(ms,options,vbcs[ii])
        dndm,dummydm= calc_dndm(ms,transfer.ks,options,vbcs[ii],silent=False)
        dndmv= dndm/dmstardm/dmvdmstar
        indx= (ms < 10.**13.)
        mv= mv[indx]
        dndmv= dndmv[indx]
        thisInterp= interpolate.InterpolatedUnivariateSpline(mv,numpy.log(dndmv))
        avg+= numpy.exp(thisInterp(mvs))*pvbc(vbcs[ii])
        print avg
        norm+= pvbc(vbcs[ii])
    return (mvs, avg/norm)

def plot_lumfunc(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.plotfile is None:
        parser.print_help()
        return
    #Calculate the halo mass function for vbc=0, 1, 2, 4 sigma
    if os.path.exists(args[0]):
        #Saved, restore
        savefile= open(args[0],'rb')
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
        raise IOError(args[0]+" does not exist; calculate using plot_massfunction.py first")
    if os.path.exists(args[1]):
        #Saved, restore
        savefile= open(args[1],'rb')
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
        #Calculate!
        #For each vbc calculate mf, mstar, mv, and jacobians
        print "Working on vbc= 0 ..."
        mstar0, mv0, dmstardm0, dmvdmstar0, mf0, rlss0= calc_mstaretc(ms,options,0.)
        print "Working on vbc= 1 ..."
        mstar1, mv1, dmstardm1, dmvdmstar1, mf1, rlss1= calc_mstaretc(ms,options,
                                                          transfer.sbc)
        print "Working on vbc= 2 ..."
        mstar2, mv2, dmstardm2, dmvdmstar2, mf2, rlss2= calc_mstaretc(ms,options,
                                                          2.*transfer.sbc)
        print "Working on vbc= 4 ..."
        mstar4, mv4, dmstardm4, dmvdmstar4, mf4, rlss4= calc_mstaretc(ms,options,
                                                          3.*transfer.sbc)
        #Save
        savefile= open(args[1],'wb')
        pickle.dump(mstar0,savefile)
        pickle.dump(mv0,savefile)
        pickle.dump(dmstardm0,savefile)
        pickle.dump(dmvdmstar0,savefile)
        pickle.dump(mf0,savefile)
        pickle.dump(rlss0,savefile)
        pickle.dump(mstar1,savefile)
        pickle.dump(mv1,savefile)
        pickle.dump(dmstardm1,savefile)
        pickle.dump(dmvdmstar1,savefile)
        pickle.dump(mf1,savefile)
        pickle.dump(rlss1,savefile)
        pickle.dump(mstar2,savefile)
        pickle.dump(mv2,savefile)
        pickle.dump(dmstardm2,savefile)
        pickle.dump(dmvdmstar2,savefile)
        pickle.dump(mf2,savefile)
        pickle.dump(rlss2,savefile)
        pickle.dump(mstar4,savefile)
        pickle.dump(mv4,savefile)
        pickle.dump(dmstardm4,savefile)
        pickle.dump(dmvdmstar4,savefile)
        pickle.dump(mf4,savefile)
        pickle.dump(rlss4,savefile)
        savefile.close()
    #Calculate dNdmv
    dndmvv0= dndmv0/dmstardm0/dmvdmstar0
    dndmvv1= dndmv1/dmstardm1/dmvdmstar1
    dndmvv2= dndmv2/dmstardm2/dmvdmstar2
    dndmvv4= dndmv4/dmstardm4/dmvdmstar4
    indx= (ms < 10.**13.)
    mv0= mv0[indx]
    mv1= mv1[indx]
    mv2= mv2[indx]
    mv4= mv4[indx]
    dndmvv0= dndmvv0[indx]
    dndmvv1= dndmvv1[indx]
    dndmvv2= dndmvv2[indx]
    dndmvv4= dndmvv4[indx]
    #Apply fs
    mv0= mv0+2.5*numpy.log10(0.01/options.fs)
    mv1= mv1+2.5*numpy.log10(0.01/options.fs)
    mv2= mv2+2.5*numpy.log10(0.01/options.fs)
    mv4= mv4+2.5*numpy.log10(0.01/options.fs)
    dndmvv0= dndmvv0*0.01/options.fs
    dndmvv1= dndmvv1*0.01/options.fs
    dndmvv2= dndmvv2*0.01/options.fs
    dndmvv4= dndmvv4*0.01/options.fs
    if options.addavg:
        if os.path.exists(args[2]):
            #Saved, restore
            savefile= open(args[2],'rb')
            mvavg= pickle.load(savefile)
            dndmvavg= pickle.load(savefile)
            savefile.close()
        else:
            #Calculate!
            mvavg, dndmvavg= calc_avg(ms,options)
            savefile= open(args[2],'wb')
            pickle.dump(mvavg,savefile)
            pickle.dump(dndmvavg,savefile)
            savefile.close()
        mvavg+= 2.5*numpy.log10(0.01/options.fs)
        dndmvavg*= 0.01/options.fs
    #Plot
    bovy_plot.bovy_print(fig_height=7.)
    #Set up top and bottom panels
    left, bottom, width, height= 0.1, 0.35, 0.8, 0.5
    axTop= pyplot.axes([left,bottom,width,height])
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.25
    axBottom= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    fig.sca(axTop)
    line1= bovy_plot.bovy_plot(mv0,dndmvv0,'k:',overplot=True,semilogy=True)
    pyplot.xlim(0.,-10.)
    axTop.set_ylim(.011,700.)
    bovy_plot._add_ticks(yticks=False)
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)    
    pyplot.ylabel(r'$\mathrm{d} N / \mathrm{d} M_V\ [\mathrm{Mpc}^{-3}]$')
    line2= bovy_plot.bovy_plot(mv1,dndmvv1,'k-',overplot=True)
    line3= bovy_plot.bovy_plot(mv2,dndmvv2,'k--',overplot=True)
    line4= bovy_plot.bovy_plot(mv4,dndmvv4,'k-.',overplot=True)
    if options.addavg:
        line5= bovy_plot.bovy_plot(mvavg,dndmvavg,
                                   '-',color='0.65',lw=2.,overplot=True)
    bovy_plot.bovy_text(r'$z=%i$' % int(options.zrei)
                        +'\n'+r'$f_s = %.2f$' % options.fs,
                        top_right=True)
    #Legend
    if options.addavg:
        pyplot.legend((line1[0],line2[0],line3[0],line4[0],line5[0]),
                      (r'$v_{bc} = 0$',
                       r'$v_{bc} = 1\,\sigma_{bc}$',
                       r'$v_{bc} = 2\,\sigma_{bc}$',
                       r'$v_{bc} = 3\,\sigma_{bc}$',
                       r'$\mathrm{averaged\ over}\ p(v_{bc})$'),
                      loc='lower left',#bbox_to_anchor=(.91,.375),
                      numpoints=2,
                      prop={'size':14},
                      frameon=False)
    else:
        pyplot.legend((line1[0],line2[0],line3[0],line4[0]),
                      (r'$v_{bc} = 0$',
                       r'$v_{bc} = 1\,\sigma_{bc}$',
                       r'$v_{bc} = 2\,\sigma_{bc}$',
                       r'$v_{bc} = 3\,\sigma_{bc}$'),
                      loc='lower left',#bbox_to_anchor=(.91,.375),
                      numpoints=2,
                      prop={'size':14},
                      frameon=False)
    #Save
    #bovy_plot.bovy_end_print(options.plotfile)              
    #return None
    #Lower panel, ratio of w/ vbc and vbc=0
    fig.sca(axBottom)
    #First interpolate to get the difference
    interpv0= interpolate.InterpolatedUnivariateSpline(mv0,numpy.log(dndmvv0))
    interpv1= interpolate.InterpolatedUnivariateSpline(mv1,numpy.log(dndmvv1))
    interpv2= interpolate.InterpolatedUnivariateSpline(mv2,numpy.log(dndmvv2))
    interpv4= interpolate.InterpolatedUnivariateSpline(mv4,numpy.log(dndmvv4))
    mvs= numpy.linspace(0.,-10.,1001)
    bovy_plot.bovy_plot([0.,-10.],[1.,1.],
                        'k:',overplot=True)
    bovy_plot.bovy_plot(mvs,numpy.exp(interpv1(mvs)-interpv0(mvs)),
                        'k-',overplot=True)
    bovy_plot.bovy_plot(mvs,numpy.exp(interpv2(mvs)-interpv0(mvs)),
                        'k--',overplot=True)
    bovy_plot.bovy_plot(mvs,numpy.exp(interpv4(mvs)-interpv0(mvs)),
                        'k-.',overplot=True)
    if options.addavg:
        print mvavg, dndmvavg
        bovy_plot.bovy_plot(mvavg,dndmvavg/numpy.exp(interpv0(mvavg)),
                            '-',color='0.65',lw=2.,overplot=True)
    #Add fs=0.01 arrow
    dmv= 2.5*numpy.log10(options.fs/0.002)
    ax=pyplot.gca()
    ax.add_patch(FancyArrowPatch((-9.,0.3),(-9.+dmv,0.3),
                                 arrowstyle='->',mutation_scale=20,fill=True,
                                 lw=1.25))
    bovy_plot.bovy_text(-7.2,.375,r'$f_s= 0.002$',
                        size=12.)
    pyplot.xlabel(r'$M_V$')
    pyplot.xlim(0.,-10.)
    axBottom.set_ylim(0.,1.199999)
    bovy_plot._add_ticks()
    pyplot.ylabel(r'$\mathrm{fractional\ effect\ wrt}\ v_{bc} = 0$')
    #Save
    bovy_plot.bovy_end_print(options.plotfile)              

def get_options():
    usage = "usage: %prog [options] <savefilename> <savefilename>\n\nsavefilename= name of the file that stuff will be saved to from plot_massfunction\nsavefilename= name of the file that lumfunc stuff will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    parser.add_option("--zrei",dest='zrei',default=15.,type='float',
                      help="Reionization redshift")
    parser.add_option("--ncost",dest='ncost',default=21,type='int',
                      help="Number of cost to use to average")
    parser.add_option("--fs",dest='fs',default=.01,type='float',
                      help="Star formation efficiency")
    parser.add_option("--addavg",action="store_true", dest="addavg",
                      default=False,
                      help="Add the average effect")
    parser.add_option("--nvbcavg",dest='nvbcavg',default=16,type='int',
                      help="Number of vbc to use to average")
    return parser

if __name__ == '__main__':
    plot_lumfunc(get_options())
