import sys
import os, os.path
from optparse import OptionParser
import cPickle as pickle
import numpy
import transfer
import massfunction
try:
    from galpy.util import bovy_plot
except ImportError:
    import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
_ERASESTR= "                                                                                "
def pvbc(vbc):
    return (3./2./numpy.pi/transfer.sbc**2.)**1.5*4.*numpy.pi*vbc**2.*numpy.exp(-1.5*vbc**2./transfer.sbc**2.)

def calc_dndm(ms,dmks,options,vbc,silent=False):
    if vbc == 0.: ncost= 1
    else: ncost= options.ncost
    #First calculate matter transfer function
    dm= numpy.zeros(len(dmks))
    for ii in range(len(dmks)):
        if not silent:
            sys.stdout.write('\r'+"Working on %i / %i ..." % (ii+1,len(dmks)))
            sys.stdout.flush()
        dm[ii]= transfer.dmAvg(options.zrei,dmks[ii],vbc,ncost=ncost)
    if not silent:
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
    #Calculate thisS
    thisS= numpy.zeros((len(ms)))
    for ii in range(len(ms)):
        thisS[ii]= massfunction.S(transfer.mToR(ms[ii]),
                                  options.zrei,vbc,ncost=ncost,
                                  window='tophat',
                                  dm=dm,
                                  dmks=dmks)
    #Calculate dndm 
    dndm= massfunction.dndm(ms,options.zrei,vbc,ncost=ncost,thisS=thisS)
    return (dndm,dm)

def calc_avg_dndm(ms,options):
    vbcs= numpy.linspace(0.,2.5,options.nvbcavg)*transfer.sbc
    avg= numpy.zeros_like(ms)
    norm= 0.
    for ii in range(1,options.nvbcavg):
        print "Working on %i / %i ..." % (ii+1,options.nvbcavg)
        dndm,dummydm= calc_dndm(ms,transfer.ks,options,vbcs[ii],silent=False)
        avg+= dndm*pvbc(vbcs[ii])
        norm+= pvbc(vbcs[ii])
    return avg/norm

def plot_massfunction(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.plotfile is None:
        parser.print_help()
        return
    #Calculate the halo mass function for vbc=0, 1, 2, 3 sigma
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
        #Calculate!
        ms= transfer.kfToMF(transfer.ks)
        dmks= transfer.ks
        dndmv0,dm0= calc_dndm(ms,dmks,options,0.)
        dndmv1,dm1= calc_dndm(ms,dmks,options,transfer.sbc)
        dndmv2,dm2= calc_dndm(ms,dmks,options,2.*transfer.sbc)
        dndmv4,dm4= calc_dndm(ms,dmks,options,3.*transfer.sbc)
        #Save
        savefile= open(args[0],'wb')
        pickle.dump(ms,savefile)
        pickle.dump(dmks,savefile)
        pickle.dump(dndmv0,savefile)
        pickle.dump(dm0,savefile)
        pickle.dump(dndmv1,savefile)
        pickle.dump(dm1,savefile)
        pickle.dump(dndmv2,savefile)
        pickle.dump(dm2,savefile)
        pickle.dump(dndmv4,savefile)
        pickle.dump(dm4,savefile)
        savefile.close()
    if options.addavg:
        if os.path.exists(args[1]):
            #Saved, restore
            savefile= open(args[1],'rb')
            dndmavg= pickle.load(savefile)
            savefile.close()
        else:
            #Calculate!
            dndmavg= calc_avg_dndm(ms,options)
            savefile= open(args[1],'wb')
            pickle.dump(dndmavg,savefile)
            savefile.close()
    #Plot
    bovy_plot.bovy_print(fig_height=7.)
    #Set up top and bottom panels
    left, bottom, width, height= 0.1, 0.35, 0.8, 0.5
    axTop= pyplot.axes([left,bottom,width,height])
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.25
    axBottom= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    fig.sca(axTop)
    line1= bovy_plot.bovy_plot(ms,dndmv0*ms,'k:',overplot=True,loglog=True)
    pyplot.xlim(10.**4.,10.**11.)
    axTop.set_ylim(10.**-6.9,10.**6.6)
    #bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)    
    pyplot.ylabel(r'$\mathrm{d} N / \mathrm{d} \ln M_{\mathrm{halo}}\ [\mathrm{Mpc}^{-3}]$')
    line2= bovy_plot.bovy_plot(ms,dndmv1*ms,'k-',overplot=True)
    line3= bovy_plot.bovy_plot(ms,dndmv2*ms,'k--',overplot=True)
    line4= bovy_plot.bovy_plot(ms,dndmv4*ms,'k-.',overplot=True)
    if options.addavg:
        line5= bovy_plot.bovy_plot(ms,dndmavg*ms,
                                   '-',color='0.65',lw=2.,overplot=True)
    bovy_plot.bovy_text(r'$z=%i$' % int(options.zrei),top_right=True)
    #Legend
    if options.addavg:
        pyplot.legend((line1[0],line2[0],line3[0],line4[0],line5[0]),
                      (r'$v_{bc} = 0$',
                       r'$v_{bc} = 1\,\sigma_{bc}, P(>v_{bc}) = 0.392$',
                       r'$v_{bc} = 2\,\sigma_{bc}, P(>v_{bc}) = 7\times10^{-3}$',
                       r'$v_{bc} = 3\,\sigma_{bc}, P(>v_{bc}) = 6\times10^{-6}$',
                       r'$\mathrm{averaged\ over}\ p(v_{bc})$'),
                      loc='lower left',#bbox_to_anchor=(.91,.375),
                      numpoints=2,
                      prop={'size':14},
                      frameon=False)
    else:
        pyplot.legend((line1[0],line2[0],line3[0],line4[0]),
                      (r'$v_{bc} = 0$',
                       r'$v_{bc} = 1\,\sigma_{bc}, P(>v_{bc}) = 0.392$',
                       r'$v_{bc} = 2\,\sigma_{bc}, P(>v_{bc}) = 7\times10^{-3}$',
                       r'$v_{bc} = 3\,\sigma_{bc}, P(>v_{bc}) = 6\times10^{-6}$'),
                      loc='lower left',#bbox_to_anchor=(.91,.375),
                      numpoints=2,
                      prop={'size':14},
                      frameon=False)
    #Lower panel, ratio of w/ vbc and vbc=0
    fig.sca(axBottom)
    bovy_plot.bovy_plot([10.**4.,10.**11],[1.,1.],
                        'k:',semilogx=True,overplot=True)
    bovy_plot.bovy_plot(ms,dndmv1/dndmv0,'k-',semilogx=True,overplot=True)
    bovy_plot.bovy_plot(ms,dndmv2/dndmv0,'k--',semilogx=True,overplot=True)
    bovy_plot.bovy_plot(ms,dndmv4/dndmv0,'k-.',semilogx=True,overplot=True)
    if options.addavg:
        bovy_plot.bovy_plot(ms,dndmavg/dndmv0,
                            '-',color='0.65',lw=2.,overplot=True)
    pyplot.xlabel(r'$M_{\mathrm{halo}} / M_\odot$')
    pyplot.xlim(10.**4.,10.**11.)
    axBottom.set_ylim(0.,1.199999)
    bovy_plot._add_ticks(xticks=False)
    pyplot.ylabel(r'$\mathrm{fractional\ effect\ wrt}\ v_{bc} = 0$')
    #Save
    bovy_plot.bovy_end_print(options.plotfile)              

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that stuff will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    parser.add_option("--zrei",dest='zrei',default=15.,type='float',
                      help="Reionization redshift")
    parser.add_option("--ncost",dest='ncost',default=21,type='int',
                      help="Number of cost to use to average")
    parser.add_option("--addavg",action="store_true", dest="addavg",
                      default=False,
                      help="Add the average effect")
    parser.add_option("--nvbcavg",dest='nvbcavg',default=16,type='int',
                      help="Number of vbc to use to average")
    return parser

if __name__ == '__main__':
    plot_massfunction(get_options())
