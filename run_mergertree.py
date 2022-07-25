#Run a mergertree simulation, for vbc=0,1,2,3 sigma using inputs from plot_massfunction
import sys
import os, os.path
from optparse import OptionParser
import cPickle as pickle
import numpy
import mergertree
def run_mergertree(parser):
    (options,args)= parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        return
    #Check whether the savefile exists
    if os.path.exists(args[0]):
        print "Savefile "+args[0]+" exists, *not* re-running ..."
    numpy.random.seed(options.seed)
    #Restore the halo mass function for vbc=0, 1, 2, 3 sigma
    if os.path.exists(args[1]):
        #Saved, restore
        savefile= open(args[1],'rb')
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
        raise IOError("mass function savefile not found, calculate first with plot_massfunction.py")
    #Load correct dm
    if options.vbc == 0:
        dm= dm0
    elif options.vbc == 1:
        dm= dm1
    elif options.vbc == 2:
        dm= dm2
    elif options.vbc == 3:
        dm= dm4
    #Run mergertree
    out= mergertree.mergertree(mres=options.mres,minit=options.minit,
                               zend=options.zrei,dzsub=0.01,
                               dm=dm,dmks=dmks)
    zs, out, zacc, newhalos, newhalosprerei, zcmain, subzc= out    
    #Save
    savefile= open(args[0],'wb')
    pickle.dump(zs,savefile)
    pickle.dump(out,savefile)
    pickle.dump(zacc,savefile)
    pickle.dump(newhalos,savefile)
    pickle.dump(newhalosprerei,savefile)
    pickle.dump(zcmain,savefile)
    pickle.dump(subzc,savefile)
    savefile.close()
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename> <savefilename>\n\nsavefilename= name of the file that stuff will be saved to\nsavefilename= savefile from plot_massfucntion.py"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("--zrei",dest='zrei',default=15.,type='float',
                      help="Reionization redshift")
    parser.add_option("--vbc",dest='vbc',default=0,type='int',
                      help="Vbc in sigma (0,1,2,3)")
    parser.add_option("--seed",dest='seed',default=1,type='int',
                      help="Seed for random number generator")
    parser.add_option("--mres",dest='mres',default=10**5.,type='float',
                      help="Mass resolution of the tree")
    parser.add_option("--minit",dest='minit',default=10**12.,type='float',
                      help="Mass of main halo")
    return parser

if __name__ == '__main__':
    run_mergertree(get_options())
