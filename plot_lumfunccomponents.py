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
_ERASESTR= "                                                                                "
def plot_lumfunccomponents(parser):
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
        raise IOError(args[1]+" does not exist; calculate using plot_lumfunc.py first")
    #Calculate dNdmv, mass function
    dndmvv0= dndmv0/dmstardm0/dmvdmstar0
    dndmvv1= dndmv1/dmstardm0/dmvdmstar0
    dndmvv2= dndmv2/dmstardm0/dmvdmstar0
    dndmvv4= dndmv4/dmstardm0/dmvdmstar0
    indx= (ms < 10.**13.)
    mv0plot= mv0[indx]
    mv1plot= mv0[indx]
    mv2plot= mv0[indx]
    mv4plot= mv0[indx]
    dndmvv0= dndmvv0[indx]
    dndmvv1= dndmvv1[indx]
    dndmvv2= dndmvv2[indx]
    dndmvv4= dndmvv4[indx]
    #Apply fs
    mv0plot+= 2.5*numpy.log10(0.01/options.fs)
    mv1plot+= 2.5*numpy.log10(0.01/options.fs)
    mv2plot+= 2.5*numpy.log10(0.01/options.fs)
    mv4plot+= 2.5*numpy.log10(0.01/options.fs)
    dndmvv0= dndmvv0*0.01/options.fs
    dndmvv1= dndmvv1*0.01/options.fs
    dndmvv2= dndmvv2*0.01/options.fs
    dndmvv4= dndmvv4*0.01/options.fs
    if options.addavg:
        if os.path.exists(args[2]):
            #Saved, restore
            savefile= open(args[2],'rb')
            dndmavg= pickle.load(savefile)
            savefile.close()
        else:
            raise IOError(args[2]+" does not exist; calculate using plot_massfunction.py first, w/ addavg")
        mvavg= mv0plot
        dndmvavg= dndmavg/dmstardm0/dmvdmstar0
        mvavg+= 2.5*numpy.log10(0.01/options.fs)
        dndmvavg*= 0.01/options.fs
        dndmvavg= dndmvavg[indx]
    #Plot
    bovy_plot.bovy_print(fig_height=4.)
    #Set up top and bottom panels
    left, bottom, width, height= 0.1, 0.5, 0.8, 0.4
    axTop= pyplot.axes([left,bottom,width,height])
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.4
    axBottom= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    fig.sca(axTop)
    bovy_plot.bovy_plot([0.,-10.],[1.,1.],'k:',
                         overplot=True)
    pyplot.xlim(0.,-10.)
    axTop.set_ylim(0.,1.199999)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)    
    bovy_plot.bovy_plot(mv1plot,dndmvv1/dndmvv0,'k-',overplot=True)
    bovy_plot.bovy_plot(mv2plot,dndmvv2/dndmvv0,'k--',overplot=True)
    bovy_plot.bovy_plot(mv4plot,dndmvv4/dndmvv0,'k-.',overplot=True)
    bovy_plot.bovy_text(r'$z=%i$' % int(options.zrei),
                        top_left=True)
    bovy_plot.bovy_text(r'$\mathrm{contribution\ from\ mass\ function}$',
                        top_right=True)
    #Save
    #bovy_plot.bovy_end_print(options.plotfile)              
    #return None
    #Lower panel, ratio of w/ vbc and vbc=0
    #Calculate dNdmv, gas fraction
    dndmvv0= dndmv0/dmstardm0/dmvdmstar0
    dndmvv1= dndmv0/dmstardm1/dmvdmstar1
    dndmvv2= dndmv0/dmstardm2/dmvdmstar2
    dndmvv4= dndmv0/dmstardm4/dmvdmstar4
    indx= (ms < 10.**13.)
    mv0plot= mv0[indx]
    mv1plot= mv1[indx]
    mv2plot= mv2[indx]
    mv4plot= mv4[indx]
    dndmvv0= dndmvv0[indx]
    dndmvv1= dndmvv1[indx]
    dndmvv2= dndmvv2[indx]
    dndmvv4= dndmvv4[indx]
    #Apply fs
    mv0plot+= 2.5*numpy.log10(0.01/options.fs)
    mv1plot+= 2.5*numpy.log10(0.01/options.fs)
    mv2plot+= 2.5*numpy.log10(0.01/options.fs)
    mv4plot+= 2.5*numpy.log10(0.01/options.fs)
    dndmvv0= dndmvv0*0.01/options.fs
    dndmvv1= dndmvv1*0.01/options.fs
    dndmvv2= dndmvv2*0.01/options.fs
    dndmvv4= dndmvv4*0.01/options.fs
    if options.addavg:
        if os.path.exists(args[3]):
            #Saved, restore
            savefile= open(args[3],'rb')
            mvavg= pickle.load(savefile)
            dndmvavg= pickle.load(savefile)
            savefile.close()
        else:
            raise IOError(args[3]+" does not exist; calculate using plot_lumfunc.py first, w/ addavg")
        #Interpolate to put on same grid
    fig.sca(axBottom)
    #First interpolate to get the difference
    interpv0= interpolate.InterpolatedUnivariateSpline(mv0plot,
                                                       numpy.log(dndmvv0))
    interpv1= interpolate.InterpolatedUnivariateSpline(mv1plot,
                                                       numpy.log(dndmvv1))
    interpv2= interpolate.InterpolatedUnivariateSpline(mv2plot,
                                                       numpy.log(dndmvv2))
    interpv4= interpolate.InterpolatedUnivariateSpline(mv4plot,
                                                       numpy.log(dndmvv4))
    mvs= numpy.linspace(0.,-10.,1001)
    line1= bovy_plot.bovy_plot([0.,-10.],[1.,1.],
                               'k:',overplot=True)
    line2= bovy_plot.bovy_plot(mvs,numpy.exp(interpv1(mvs)-interpv0(mvs)),
                               'k-',overplot=True)
    line3= bovy_plot.bovy_plot(mvs,numpy.exp(interpv2(mvs)-interpv0(mvs)),
                               'k--',overplot=True)
    line4= bovy_plot.bovy_plot(mvs,numpy.exp(interpv4(mvs)-interpv0(mvs)),
                               'k-.',overplot=True)
    if options.addavg:
        line5= bovy_plot.bovy_plot(mvavg,dndmvavg/numpy.exp(interpv0(mvavg)),
                                   '-',color='0.65',lw=2.,overplot=True)
    bovy_plot.bovy_text(r'$\mathrm{contribution\ from\ gas\ fraction}$',
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
                      loc='lower right',#bbox_to_anchor=(.91,.375),
                      numpoints=2,
                      prop={'size':14},
                      frameon=False)
    pyplot.xlabel(r'$M_V$')
    pyplot.xlim(0.,-10.)
    axBottom.set_ylim(0.,1.2499999)
    bovy_plot._add_ticks()
    pyplot.text(1.25,2.25,r'$\mathrm{fractional\ effect\ in}\ \mathrm{d} N / \mathrm{d} M_V\ \mathrm{wrt}\ v_{bc} = 0$',rotation=90.,size=14.)
    #pyplot.ylabel(r'$\mathrm{fractional\ diff\ wrt}\ v_{bc} = 0$')
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
    plot_lumfunccomponents(get_options())
