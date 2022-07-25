import numpy
class vialactea:
    """Holds data from VL"""
    def __init__(self):
        #Read everything
        self.GCdistance= numpy.loadtxt('../vialactea/progGCdistance.txt')
        self.Mtidal= numpy.loadtxt('../vialactea/progMtidal.txt')
        self.rtidal= numpy.loadtxt('../vialactea/progrtidal.txt')
        self.Vmax= numpy.loadtxt('../vialactea/progVmax.txt')
        self.rVmax= numpy.loadtxt('../vialactea/progrVmax.txt')
        self.X= numpy.loadtxt('../vialactea/progX.txt')*40. #convert to coMpc
        self.Y= numpy.loadtxt('../vialactea/progY.txt')*40.
        self.Z= numpy.loadtxt('../vialactea/progZ.txt')*40.
        self.VX= numpy.loadtxt('../vialactea/progVX.txt')
        self.VY= numpy.loadtxt('../vialactea/progVY.txt')
        self.VZ= numpy.loadtxt('../vialactea/progVZ.txt')
        #Read z=0 snapshot
        self._zzero= numpy.loadtxt('../vialactea/vltwosubs.txt')
        self.GCdistance[:,-1]= self._zzero[:,1]
        self.Mtidal[:,-1]= self._zzero[:,5]
        self.rtidal[:,-1]= self._zzero[:,6]
        self.Vmax[:,-1]= self._zzero[:,3]
        self.rVmax[:,-1]= self._zzero[:,4]
        self.X[:,-1]= self._zzero[:,7]/1000. #convert to Mpc
        self.Y[:,-1]= self._zzero[:,8]/1000.
        self.Z[:,-1]= self._zzero[:,9]/1000.
        self.VX[:,-1]= self._zzero[:,10]
        self.VY[:,-1]= self._zzero[:,11]
        self.VZ[:,-1]= self._zzero[:,12]
        #Center X,Y,Z on host
        for ii in range(len(self.X[0,:])):
            self.X[:,ii]= self.X[:,ii]-self.X[-2,ii]
            self.Y[:,ii]= self.Y[:,ii]-self.Y[-2,ii]
            self.Z[:,ii]= self.Z[:,ii]-self.Z[-2,ii]
        #Read redshifts
        self.redshifts= numpy.loadtxt('../vialactea/stepToTimeVL2.txt')[:,2]
        return None

    def mv(self,fs=0.05,fgas=0.17,mvpersolarmass=6.7):
        out= 6.7-2.5*numpy.log10(fs*fgas*self.Mtidal)
        out[numpy.isinf(out)]= numpy.nan
        return out

    def cumulmass(self,indx,r200=402.):
        """Calculate the cumulative mass function at the redshift corresponding to indx"""
        ms= self.Mtidal[:,indx]
        ms= ms[(self.GCdistance[:,-1] < r200)]
        ms= sorted(ms)[::-1]
        cumulmass= numpy.ones_like(ms)
        return (ms,numpy.cumsum(cumulmass))
