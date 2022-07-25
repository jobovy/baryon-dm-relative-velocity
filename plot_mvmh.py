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
def plot_mvmh(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.plotfile is None:
        parser.print_help()
        return
    #Restore
    if os.path.exists(args[0]):
        #Saved, restore
        savefile= open(args[0],'rb')
        ms= pickle.load(savefile)
        dndmv0= pickle.load(savefile)
        dndmv1= pickle.load(savefile)
        dndmv2= pickle.load(savefile)
        dndmv4= pickle.load(savefile)
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
        mstar1= pickle.load(savefile)
        mv1= pickle.load(savefile)
        dmstardm1= pickle.load(savefile)
        dmvdmstar1= pickle.load(savefile)
        mstar2= pickle.load(savefile)
        mv2= pickle.load(savefile)
        dmstardm2= pickle.load(savefile)
        dmvdmstar2= pickle.load(savefile)
        mstar4= pickle.load(savefile)
        mv4= pickle.load(savefile)
        dmstardm4= pickle.load(savefile)
        dmvdmstar4= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError(args[1]+" does not exist; calculate using plot_lumfunc.py first")
    bovy_plot.bovy_print()
    line1= bovy_plot.bovy_plot(mv0,ms,'k:',overplot=False,semilogy=True,
                               xlabel=r'$M_V\ [\mathrm{mag}]$',
                               ylabel=r'$M_\mathrm{halo} / M_\odot$',
                               xrange=[0.,-10.],
                               yrange=[10.**5,10.**10.])
    line2= bovy_plot.bovy_plot(mv1,ms,'k-',overplot=True,semilogy=True)
    line3= bovy_plot.bovy_plot(mv2,ms,'k--',overplot=True,semilogy=True)
    line4= bovy_plot.bovy_plot(mv4,ms,'k-.',overplot=True,semilogy=True)
    bovy_plot.bovy_text(r'$z=%i$' % int(options.zrei)
                        +'\n'+r'$f_s = %.3f$' % options.fs,
                        top_left=True)
    #Legend
    pyplot.legend((line1[0],line2[0],line3[0],line4[0]),
                  (r'$v_{bc} = 0$',
                   r'$v_{bc} = 1\,\sigma_{bc}$',
                   r'$v_{bc} = 2\,\sigma_{bc}$',
                   r'$v_{bc} = 4\,\sigma_{bc}$'),
                  loc='lower right',#bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':14},
                  frameon=False)
    #Save
    bovy_plot.bovy_end_print(options.plotfile)              
    return None
    
    
def get_options():
    usage = "usage: %prog [options] <savefilename> <savefilename>\n\nsavefilename= name of the file that stuff will be saved to from plot_massfunction\nsavefilename= name of the file that lumfunc stuff will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    parser.add_option("--zrei",dest='zrei',default=15.,type='float',
                      help="Reionization redshift")
    parser.add_option("--fs",dest='fs',default=.01,type='float',
                      help="Star formation efficiency")
    return parser

if __name__ == '__main__':
    plot_mvmh(get_options())
