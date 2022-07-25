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
def find_nearest(array,value):
    idx=(numpy.abs(array-value)).argmin()
    return idx
def plot_sdsslumfunc(parser):
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
    dev1= pickle.load(savefile)
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
                                          s=0.08*len(xs0)/4.)
    interp1= interpolate.UnivariateSpline(xs1,numpy.log(lumfunc1),
                                          s=0.08*len(xs0)/4.)
    interp2= interpolate.UnivariateSpline(xs2,numpy.log(lumfunc2),
                                          s=0.08*len(xs0)/4.)
    interp3= interpolate.UnivariateSpline(xs3,numpy.log(lumfunc3),
                                          s=0.08*len(xs0)/4.)
    mvs= numpy.linspace(-10.,0.,1001)
    loglumfunc0= interp0(mvs)
    loglumfunc1= interp1(mvs)
    loglumfunc2= interp2(mvs)
    loglumfunc3= interp3(mvs)
    #Plot
    bovy_plot.bovy_print()
    line1= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc0),
                               'k:',semilogy=True,
                               xrange=[0.,-10.],
                               yrange=[0.03,20.],
                               xlabel=r'$M_V$',
                               ylabel=r'$\mathrm{d} N / \mathrm{d} M_V$')
    line2= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc1),'k-',overplot=True)
    line3= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc2),'k--',overplot=True)
    line4= bovy_plot.bovy_plot(mvs,numpy.exp(loglumfunc3),'k-.',overplot=True)
    bovy_plot._add_ticks(yticks=False)
    #Add errorbars for lumfunc1
    errorxs= [-6.,-4.,-2.,-0.2]
    errorindxs= []
    for ii in range(len(errorxs)): errorindxs.append(find_nearest(xs1,errorxs[ii]))
    print xs1[errorindxs],dev1[:,errorindxs]
    for ii in range(len(xs1)): 
        dev1[:,ii]-= lumfunc1[ii]
        dev1[0,ii]*= -1
    print xs1[errorindxs],dev1[:,errorindxs]
    pyplot.errorbar(xs1[errorindxs],lumfunc1[errorindxs],
                    yerr=dev1[:,errorindxs],
                    fmt=',',ecolor='k',marker='None')
    #Add Koposov points
    koposov= numpy.loadtxt('koposov09a.txt',delimiter='|')
    bovy_plot.bovy_plot(koposov[:,0],koposov[:,1],'Dw',overplot=True,
                        mew=1.)
    #Upper limit
    #ax= pyplot.gca()
    #ax.add_patch(FancyArrowPatch((-0.2,0.1),(-0.2,0.05),
    #                             arrowstyle='->',mutation_scale=20,fill=True,
    #                             lw=1.25))
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
    bovy_plot.bovy_text(r'$z=0,\ 10^{12}\ M_\odot\ \mathrm{parent\ halo}$'
                        +'\n'+r'$\mathrm{SDSS\ volume\ selection}$',
                        top_right=True)
    #Save
    bovy_plot.bovy_end_print(options.plotfile)              
    return None

def get_options():
    usage = "usage: %prog [options] <savefilenames> \n\nsavefilenames= name of the files that hold the averages computed with plot_mergertree.py"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    return parser

if __name__ == '__main__':
    plot_sdsslumfunc(get_options())
