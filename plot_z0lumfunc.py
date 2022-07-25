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
from matplotlib.ticker import NullFormatter, FuncFormatter, MultipleLocator
from matplotlib.patches import FancyArrowPatch
from plot_massfunction import calc_dndm
_ERASESTR= "                                                                                "
def plot_z0lumfunc(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.plotfile is None:
        parser.print_help()
        return
    #Load all averages
    #vbc=0
    savefile= open(args[0],'rb')
    xs0= pickle.load(savefile)
    lumfunc0= pickle.load(savefile)
    savefile.close()
    #vbc=1
    savefile= open(args[1],'rb')
    xs1= pickle.load(savefile)
    lumfunc1= pickle.load(savefile)
    dev= pickle.load(savefile)
    allhists= pickle.load(savefile)
    mvmassmvs= pickle.load(savefile)
    avgmvmass= pickle.load(savefile)
    savefile.close()
    #vbc=2
    savefile= open(args[2],'rb')
    xs2= pickle.load(savefile)
    lumfunc2= pickle.load(savefile)
    savefile.close()
    #vbc=3
    savefile= open(args[3],'rb')
    xs3= pickle.load(savefile)
    lumfunc3= pickle.load(savefile)
    savefile.close()
    #Interpolate all
    interp0= interpolate.UnivariateSpline(xs0,numpy.log(lumfunc0),
                                          s=0.02*len(xs0)/4.)
    interp1= interpolate.UnivariateSpline(xs1,numpy.log(lumfunc1),
                                          s=0.02*len(xs0)/4.)
    interp2= interpolate.UnivariateSpline(xs1,numpy.log(lumfunc2),
                                          s=0.02*len(xs0)/4.)
    interp3= interpolate.UnivariateSpline(xs3,numpy.log(lumfunc3),
                                          s=0.02*len(xs0)/4.)
    #Interpolate mv vs. mass
    interpmvmass= interpolate.UnivariateSpline(((mvmassmvs+numpy.roll(mvmassmvs,-1))/2.)[0:-1],
                                               numpy.log(avgmvmass),
                                               s=0.02*len(xs0)/4.)
    mvs= numpy.linspace(-10.,0.,1001)
    loglumfunc0= interp0(mvs)
    loglumfunc1= interp1(mvs)
    loglumfunc2= interp2(mvs)
    loglumfunc3= interp3(mvs)
    #Plot
    bovy_plot.bovy_print(fig_height=7.)
    #Set up top and bottom panels
    left, bottom, width, height= 0.1, 0.35, 0.8, 0.5
    axTop= pyplot.axes([left,bottom,width,height])
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.25
    axBottom= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    fig.sca(axTop)
    line1= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc0),
                               'k:',overplot=True,semilogy=True)
    pyplot.xlim(0.,-10.)
    axTop.set_ylim(1.1,10000.)
    bovy_plot._add_ticks(yticks=False)
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)    
    pyplot.ylabel(r'$\mathrm{d} N / \mathrm{d} M_V$')
    line2= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc1),'k-',overplot=True)
    line3= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc2),'k--',overplot=True)
    line4= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc3),'k-.',overplot=True)
    #if options.addavg:
    #    line5= bovy_plot.bovy_plot(mvavg,dndmvavg,
    #                               '-',color='0.65',lw=2.,overplot=True)
    #Legend
    if False:#options.addavg:
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
    bovy_plot.bovy_text(r'$z=0,\ 10^{12}\ M_\odot\ \mathrm{parent\ halo}$',
                        top_right=True)

    #Add twin x axis
    ax= axTop
    def my_formatter(x, pos):
        mh= numpy.exp(interpmvmass(x))
        expon= numpy.log10(mh)
        if numpy.fabs(expon-round(expon)) >= 0.05:
            return r'$%i \times 10^{%.0f}$' % (10.**(expon-numpy.floor(expon)),
                                                     numpy.floor(expon))
        else:
            return r'$10^{%.0f}$' % round(expon)
    ax2= pyplot.twiny()
#    ax2.xaxis.set_major_locator(MultipleLocator(1.5))
    major_formatter = FuncFormatter(my_formatter)
    ax2.xaxis.set_major_formatter(major_formatter)
    #xstep= ax.xaxis.get_majorticklocs()
    #xstep= xstep[1]-xstep[0]
    #ax2.xaxis.set_minor_locator(MultipleLocator(xstep/5.))
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    xmin, xmax= ax.xaxis.get_view_interval()
    ax2.xaxis.set_view_interval(xmin,xmax,ignore=True)
    ax2.set_xlabel('$M_{\mathrm{halo}} / M_\odot$')

    #Save
    #bovy_plot.bovy_end_print(options.plotfile)              
    #return None
    #Lower panel, ratio of w/ vbc and vbc=0
    fig.sca(axBottom)
    bovy_plot.bovy_plot([0.,-10.],[1.,1.],
                        'k:',overplot=True)
    interp1= interpolate.UnivariateSpline(xs1,numpy.log(lumfunc1/lumfunc0),
                                          s=0.2*len(xs0)/4.)
    interp2= interpolate.UnivariateSpline(xs1,numpy.log(lumfunc2/lumfunc0),
                                          s=0.2*len(xs0)/4.)
    interp3= interpolate.UnivariateSpline(xs3,numpy.log(lumfunc3/lumfunc0),
                                          s=0.2*len(xs0)/4.)
    #bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc1-loglumfunc0),
    bovy_plot.bovy_plot(mvs,numpy.exp(interp1(mvs)),
                        'k-',overplot=True)
    #bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc2-loglumfunc0),
    bovy_plot.bovy_plot(mvs,numpy.exp(interp2(mvs)),
                        'k--',overplot=True)
    #bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc3-loglumfunc0),
    bovy_plot.bovy_plot(mvs,numpy.exp(interp3(mvs)),
                        'k-.',overplot=True)
    if False:#options.addavg:
        print mvavg, dndmvavg
        bovy_plot.bovy_plot(mvavg,dndmvavg/numpy.exp(interpv0(mvavg)),
                            '-',color='0.65',lw=2.,overplot=True)
    pyplot.xlabel(r'$M_V$')
    pyplot.xlim(0.,-10.)
    axBottom.set_ylim(0.,1.199999)
    bovy_plot._add_ticks()
    pyplot.ylabel(r'$\mathrm{fractional\ effect\ wrt}\ v_{bc} = 0$')
    #Save
    bovy_plot.bovy_end_print(options.plotfile)              
    

def get_options():
    usage = "usage: %prog [options] <savefilenames> \n\nsavefilenames= name of the files that hold the averages computed with plot_mergertree.py"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    return parser

if __name__ == '__main__':
    plot_z0lumfunc(get_options())
