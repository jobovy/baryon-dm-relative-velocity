import sys
import os, os.path
from optparse import OptionParser
import numpy
from scipy import interpolate, integrate, special
import transfer
try:
    from galpy.util import bovy_plot
except ImportError:
    import bovy_plot
def plot_correlation(parser):
    (options,args)= parser.parse_args()
    if options.plotfile is None:
        parser.print_help()
        return
    if not options.zend == 1020.:
        raise NotImplementedError("z < 1020 not implemented")
    else:
        ks= transfer.ks
        vbc2= transfer.matterHz(options.zend)**2.\
            *transfer._DELTA2ZETA*(ks/0.002)**transfer._NSMINUSONE\
            *(transfer.dDeltab-transfer.dDeltac)**2./ks**3.
        indx= (ks > 0.001)*(ks < 10.)
        ks= ks[indx]
        vbc2= vbc2[indx]
    #Interpolate function
    interpVbc2= interpolate.interp1d(ks,vbc2,kind=3)
    if True:
        bovy_plot.bovy_print()
        thisks= 10.**(numpy.linspace(numpy.log10(numpy.amin(ks))+0.01,
                                     numpy.log10(numpy.amax(ks))-0.01,1001))
        bovy_plot.bovy_plot(thisks,interpVbc2(thisks),'k-',lw=2.,loglog=True)
        bovy_plot.bovy_plot(ks,vbc2,'ko',overplot=True)
        bovy_plot.bovy_end_print('testinterp.png')
    #Integrate
    rs= 10.**(numpy.linspace(-1.,2.5,1001))
    nrs= len(rs)
    xivbc= numpy.zeros_like(rs)
    xivbck2= numpy.zeros_like(rs)
    for ii in range(nrs):
        xivbc[ii]= integrate.quad((lambda x: interpVbc2(x)\
                                       *special.sph_jn(0.,x*rs[ii])[0]),
                                  0.005,1.)[0]
        if False:
            xivbck2[ii]= integrate.quad((lambda x: x**2.*interpVbc2(x)\
                                             *special.sph_jn(0.,x*rs[ii])[0]),
                                        0.005,1.)[0]
    if False:
        print numpy.sqrt(numpy.amax(xivbc))*transfer._C
        print numpy.sqrt(numpy.amax(xivbck2))*transfer._C
    xivbc/= numpy.amax(xivbc)
    xivbck2/= numpy.amax(xivbck2)
    #Plot
    bovy_plot.bovy_print()
    if not options.loglog:
        bovy_plot.bovy_plot(rs,xivbc,'k-',lw=1.5,
                            xrange=[0.1,10.**2.5],
                            yrange=[-.2,1.],
                            loglog=False,
                            semilogx=True,
                            xlabel=r'$r\ [\mathrm{Mpc}]$',
                            ylabel=r'$\xi_{v_{bc}}(r)$',zorder=5)
        if False:
            bovy_plot.bovy_plot(rs,xivbck2,'k--',lw=1.5,overplot=True)
        bovy_plot.bovy_plot([0.1,10.**2.5],[0.,0.],'-',color='0.65',
                            overplot=True,lw=2.,zorder=1)
        bovy_plot.bovy_plot([2.,2.],[-1.,1.],'-',color='0.65',
                            overplot=True,lw=2.,zorder=2)
        bovy_plot.bovy_text(1.52,0.3,r'$\mathrm{Local\ Group}$',rotation='vertical')
        bovy_plot.bovy_plot([80./transfer._H0NOC*100.,
                             80./transfer._H0NOC*100.],
                            [-1.,1.],'-',color='0.65',
                            overplot=True,lw=2.,zorder=2)
        bovy_plot.bovy_text(80./transfer._H0NOC*100.-28.,
                            0.9,r'$\mathrm{MARK\ III}$',
                            rotation='vertical')
    else:
        posindx= (xivbc >= 0.)
        bovy_plot.bovy_plot(rs[posindx],xivbc[posindx],'k-',lw=1.5,
                            xrange=[0.1,10.**2.5],
                            yrange=[0.01,1.],
                            loglog=True,
                            xlabel=r'$r\ [\mathrm{Mpc}]$',
                            ylabel=r'$\xi_{v_{bc}}(r)$')
        negindx= (xivbc < 0.)
        bovy_plot.bovy_plot(rs[negindx],-xivbc[negindx],'k--',lw=1.5,
                            overplot=True)
    bovy_plot.bovy_end_print(options.plotfile)
                        

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("-o",dest='plotfile',default=None,
                      help="name of the file that will hold the figure")
    parser.add_option("--zend",dest='zend',default=1020.,type='float',
                      help="End redshift (<= 1020)")
    parser.add_option("--loglog",action="store_true", dest="loglog",
                      default=False,
                      help="Plot on log log scale")
    return parser

if __name__ == '__main__':
    plot_correlation(get_options())
