import os
from pyrap.tables import table
import numpy as np
import pylab
from numpy import ma
import sys
import warnings
import time
import math
import argparse

class CovWeights:
    def __init__(self,MSName,ntsol=1,SaveDataProducts=0,uvcut=[0,2000],gainfile=None,phaseonly=False,nfreqsol=1,norm=False):
        if MSName[-1]=="/":
            self.MSName=MSName[0:-1]
        else:
            self.MSName=MSName
        self.MaxCorrTime=0
        self.SaveDataProducts=SaveDataProducts
        self.ntSol=ntsol
        self.nfreqsol=nfreqsol
        self.uvcut=uvcut
        self.gainfile=gainfile
        self.phaseonly=phaseonly
        self.normalise=norm

        
    def FindWeights(self,tcorr=0,colname=""):
        ms=table(self.MSName,readonly=False)
        # open antennas
        ants=table(ms.getkeyword("ANTENNA"))
        # open antenna tables
        antnames=ants.getcol("NAME")
        nAnt=len(antnames)
        # open uvw to perform uvcut - TEST
        u,v,_=ms.getcol("UVW").T
        # load ant indices
        A0=ms.getcol("ANTENNA1")
        A1=ms.getcol("ANTENNA2")
        tarray=ms.getcol("TIME")
        nbl=np.where(tarray==tarray[0])[0].size
        warnings.filterwarnings("ignore")
        warnings.filterwarnings("default")
        if "RESIDUAL_DATA" not in ms.colnames():
            print "RESIDUAL_DATA column not in measurement set: creating it and filling with RESIDUAL_DATA"
            desc=ms.getcoldesc("CORRECTED_DATA")
            desc["name"]="RESIDUAL_DATA"
            desc['comment']=desc['comment'].replace(" ","_")
            ms.addcols(desc)
            ms.putcol("RESIDUAL_DATA","CORRECTED_DATA")
        flags=ms.getcol("FLAG")
        print "Please ensure that RESIDUAL_DATA or CORRECTED_DATA contains residual visibilities from complete skymodel subtraction."
        residualdata=ms.getcol("RESIDUAL_DATA")
        flags=ms.getcol("FLAG")
        # apply uvcut
        uvlen=np.sqrt(u**2+v**2)
        flags[uvlen>self.uvcut[1]]=1
        flags[uvlen<self.uvcut[0]]=1
        # apply flags to data
        residualdata[flags==1]=0
        # exit files gracefully
        ants.close()
        # initialise
        nChan=residualdata.shape[1]
        nPola=residualdata.shape[2]
        nt=residualdata.shape[0]/nbl
        # reshape antennas and data columns
        residualdata=residualdata.reshape((nt,nbl,nChan,nPola))
        # average residual data within calibration cells
        ### TODO: the averaging should be done later; try using a filter instead of this ###
#        times=np.array(list(sorted(set(tarray))))
        A0=A0.reshape((nt,nbl))
        A1=A1.reshape((nt,nbl))
#        ant1=np.arange(nAnt)
#        CoeffArray=np.zeros((nt,nAnt))
#        print "Begin calculating antenna-based coefficients"
#        for i,t_i in enumerate(times):
#            indexmax=min(len(times)-1,i+self.ntSol)
#            indexmin=max(0,i-self.ntSol)
#            tmin=times[indexmin]
#            tmax=times[indexmax]
#            tmask=(tarray>tmin)*(tarray<tmax)
#            for ant in ant1:
#                # set of vis for baselines ant-ant_i
#                set1=(A0==ant)
#                # set of vis for baselines ant_i-ant
#                set2=(A1==ant)
#                resmask=tmask*(set1+set2)
#                rarray=residualdata[resmask]
#                CoeffArray[i,ant] = np.mean(np.abs(rarray*rarray.conj()))
#            PrintProgress(i,nt)


        if self.ntSol>1:
            tspill=nt%self.ntSol
            nt1=nt+self.ntSol-tspill
            for i in range(nt1/self.ntSol):
                for j in range(self.nfreqsol):
                    residualdata[i*self.ntSol:(i+1)*self.ntSol,:,self.nfreqsol*j:(j+1)*nfreqsol,:]=np.mean(residualdata[i*self.ntSol:(i+1)*self.ntSol,:,:,:],axis=0)
        A0=A0.reshape((nt,nbl))
        A1=A1.reshape((nt,nbl))
        ant1=np.arange(nAnt)
        # make rms array
        darray=ms.getcol("RESIDUAL_DATA").reshape((nt,nbl,nChan,nPola))
        ms.close()
        rmsarray=np.zeros((nt,nbl,nChan,2),dtype=np.complex64)
        residuals=np.zeros_like(rmsarray,dtype=np.complex64)
        rmsarray[:,:,:,0]=darray[:,:,:,1]
        rmsarray[:,:,:,1]=darray[:,:,:,2]
        # make proper residual array
        residuals[:,:,:,0]=darray[:,:,:,0]
        residuals[:,:,:,1]=darray[:,:,:,3]
        # antenna coefficient array
        CoeffArray=np.zeros((nt,nAnt,2))
        # start calculating the weights
        print "Begin calculating antenna-based coefficients"
        warnings.filterwarnings("ignore")
        print "Find variance-only weights"
        for t_i in range(nt):
            # build weights for each antenna at time t_i
            for ant in ant1:
                # set of vis for baselines ant-ant_i
                set1=np.where(A0[t_i]==ant1)[0]
                # set of vis for baselines ant_i-ant
                set2=np.where(A1[t_i]==ant)[0]
                CoeffArray[t_i,ant,0] = (np.mean(np.append(residuals[t_i,set1,:,:],residuals[t_i,set2,:,:])*np.append(residuals[t_i,set1,:,:],residuals[t_i,set2,:,:]).conj())\
                                              - np.std( (np.append(rmsarray[t_i,set1,:,:], rmsarray[t_i,set2,:,:]))) )
                CoeffArray[t_i,ant,1] = np.mean(np.append(residuals[t_i,set1,:,:],residuals[t_i,set2,:,:]))*(np.mean(np.append(residuals[t_i,set1,:,:],residuals[t_i,set2,:,:]))).conj()
            PrintProgress(t_i,nt)
        warnings.filterwarnings("default")
        for i in range(nAnt):
            # get rid of NaN
            CoeffArray[np.isnan(CoeffArray)]=np.inf
            tempars=CoeffArray[:,i,0]
            gtempars=CoeffArray[:,i,1]
            thres=0.25*np.median(tempars)
            gthres=0.25*np.median(gtempars)
            CoeffArray[:,i,0][tempars<thres]=thres
            CoeffArray[:,i,0][gtempars<gthres]=gthres
            # normalise per antenna
            if self.normalise==True:
                CoeffArray[:,i,0]=CoeffArray[:,i,0]/CoeffArray[:,i,1]**2
#            CoeffArray[:,i,1]=CoeffArray[:,i,1]/np.mean(CoeffArray[:,i,1])
        if colname=="":
            coeffFilename=self.MSName+"/CoeffArray.ntsol%i.npy"%(ntsol)
        else:
            coeffFilename=self.MSName+"/CoeffArray.%s.ntsol%i.npy"%(colname,ntsol)
        print "Save coefficient array as %s."%coeffFilename
        np.save(coeffFilename,CoeffArray)
        return CoeffArray
                        
    def SaveWeights(self,CoeffArray,colname=None,AverageOverChannels=True,tcorr=0):
        print "Begin saving the data"
        ms=table(self.MSName,readonly=False)
        # open antennas
        ants=table(ms.getkeyword("ANTENNA"))
        # open antenna tables
        antnames=ants.getcol("NAME")
        nAnt=len(antnames)
        tarray=ms.getcol("TIME")
        darray=ms.getcol("DATA")
        tvalues=np.array(sorted(list(set(tarray))))
        nt=tvalues.shape[0]
        nbl=tarray.shape[0]/nt
        nchan=darray.shape[1]
        A0=np.array(ms.getcol("ANTENNA1").reshape((nt,nbl)))
        A1=np.array(ms.getcol("ANTENNA2").reshape((nt,nbl)))
        if colname in ms.colnames():
            print "%s column already present; will overwrite"%colname
        else:
            W=np.ones((nt*nbl,nchan))
            desc=ms.getcoldesc("IMAGING_WEIGHT")
            desc["name"]=colname
            desc['comment']=desc['comment'].replace(" ","_")
            ms.addcols(desc)
            ms.putcol(colname,W)
        # create weight array
        w=np.zeros((nt,nbl,nchan))
        ant1=np.arange(nAnt)
        print "Fill weights array"
        A0ind=A0[0,:]
        A1ind=A1[0,:]
        warnings.filterwarnings("ignore")

        # do gains stuff
        ant1gainarray,ant2gainarray=readGainFile(self.gainfile, ms, nt, nchan, nbl,tarray,nAnt,self.MSName,self.phaseonly)
        for i in range(nbl):
            for j in range(nchan):               
                #w[:,i,j]=1./(CoeffArray[:,A0ind[i]]*ant2gainarray[i,j]+CoeffArray[:,A1ind[i]]*ant1gainarray[i,j]+CoeffArray[:,A0ind[i]]*CoeffArray[:,A1ind[i]] + 0.1)
                if self.normalise==True:
                    ga=CoeffArray[:,A0ind[i],1]**2
                    gb=CoeffArray[:,A1ind[i],1]**2
                    w[:,i,j]=1./( ga*gb*(CoeffArray[:,A0ind[i],0]+CoeffArray[:,A1ind[i],0]+CoeffArray[:,A0ind[i],0]*CoeffArray[:,A1ind[i],0]) + 0.1)
                else:
                    w[:,i,j]=1./(CoeffArray[:,A0ind[i],0]*CoeffArray[:,A1ind[i],1]+CoeffArray[:,A1ind[i],0]*CoeffArray[:,A0ind[i],1]+CoeffArray[:,A0ind[i],0]*CoeffArray[:,A1ind[i],0] + 0.1)
            PrintProgress(i,nbl)
        warnings.filterwarnings("default")
        w=w.reshape(nt*nbl,nchan)
        w[np.isnan(w)]=0
        w[np.isinf(w)]=0
        # normalise
        w=w/np.mean(w)
        # save in weights column
        if colname!=None:
            ms.putcol(colname,w)
        else: print "No colname given, so weights not saved in MS."
        ants.close()
        ms.close()

def readGainFile(gainfile,ms,nt,nchan,nbl,tarray,nAnt,msname,phaseonly):
    if phaseonly==True:
        print "Assume amplitude gain values of 1 everywhere"
        ant1gainarray1=np.ones((nt*nbl,nchan))
        ant2gainarray1=np.ones((nt*nbl,nchan))
    else:
        if gainfile[-4:]==".npz":
            print "Assume reading a kMS sols file"
            gainsnpz=np.load(gainfile)
	    gains=gainsnpz["Sols"]
            ant1gainarray=np.ones((nt*nbl,nchan))
            ant2gainarray=np.ones((nt*nbl,nchan))
            A0arr=ms.getcol("ANTENNA1")
            A1arr=ms.getcol("ANTENNA2")
            print "Build squared gain array"
            for i in range(len(gains)):
                timemask=(tarray>gains[i][0])*(tarray<gains[i][1])
                for j in range(nAnt):
                    mask1=timemask*(A0arr==j)
                    mask2=timemask*(A1arr==j)
                    for k in range(nchan):
                        ant1gainarray[mask1,:]=np.abs(np.nanmean(gains[i][2][0,j,0]))#np.abs(np.nanmean(gains[i][3][0,j]))
                        ant2gainarray[mask2,:]=np.abs(np.nanmean(gains[i][2][0,j,0]))#np.abs(np.nanmean(gains[i][3][0,j]))
                PrintProgress(i,len(gains))
            np.save(msname+"/ant1gainarray",ant1gainarray)
            np.save(msname+"/ant2gainarray",ant2gainarray)
            ant1gainarray=np.load(msname+"/ant1gainarray.npy")
            ant2gainarray=np.load(msname+"/ant2gainarray.npy")
            #        ant1gainarray1=np.ones((nt,nbl,nchan))
            #        ant2gainarray1=np.ones((nt,nbl,nchan))
            #        for i in range(nchan):
            #            ant1gainarray1[:,:,i]=ant1gainarray**2
            #            ant2gainarray1[:,:,i]=ant2gainarray**2
            ant1gainarray1=ant1gainarray**2#1.reshape((nt*nbl,nchan))
            ant2gainarray1=ant2gainarray**2#1.reshape((nt*nbl,nchan))
            if gainfile[-3:]==".h5":
                print "Assume reading losoto h5parm file"
                import losoto
                solsetName="sol000"
                soltabName="amp000"
                try:
                    gfile=losoto.h5parm.openSoltab(gainfile,solsetName=solsetName,soltabName=soltabName)
                except:
                    print "Could not find amplitude gains in h5parm. Assuming gains of 1 everywhere."
                    ant1gainarray1=np.ones((nt*nbl,nchan))
                    ant2gainarray1=np.ones((nt*nbl,nchan))
                    return ant1gainarray1,ant2gainarray1
                freqs=table(msname+"/SPECTRAL_WINDOW").getcol("CHAN_FREQ")
                gains=gfile.getValues()[0] # axes: pol, dir, ant, freq, times
                gfreqs=gfile.getValues()[1]["freq"]
                times=fgile.getValues()[1]["time"]
                ant1gainarray=np.zeros((nt*nbl,nchan))
                ant2gainarray=np.zeros((nt*nbl,nchan))
                for j in range(nAnt):
                    mask1=timemask*(A0arr==j)
                    mask2=timemask*(A1arr==j)
                    for k in range(nchan):
                        if freqs[k] in gfreqs:
                            freqmask=(gfreqs==k)
                            ant1gainarray[mask1,k]=np.mean(gains[:,0,j,freqmask],axis=0)
                            ant2gainarray[mask2,k]=np.mean(gains[:,0,j,freqmask],axis=0)
            else:
                print "Gain file type not currently supported. Assume all gain amplitudes are 1."
                ant1gainarray=np.ones((nt*nbl,nchan))
                ant2gainarray=np.ones((nt*nbl,nchan))
                                

    return ant1gainarray1,ant2gainarray1

        
### auxiliary functions ###
def PrintProgress(currentIter,maxIter,msg=""):
    sys.stdout.flush()
    if msg=="":
        msg="Progress:"
    sys.stdout.write("\r%s %5.1f %% "%(msg,100*(currentIter+1.)/maxIter))
    if currentIter==(maxIter-1):
        sys.stdout.write("\n")
def invSVD(A):
    u,s,v=np.linalg.svd(A)
    s[s<1.e-6*s.max()]=1.e-6*s.max()
    ssq=np.abs((1./s))
    # rebuild matrix
    Asq=np.dot(v,np.dot(np.diag(ssq),np.conj(u)))
    v0=v.T*ssq.reshape(1,ssq.size)
    return Asq
def readArguments():
    parser=argparse.ArgumentParser("Calculate visibility imagin weights based on calibration quality")
    parser.add_argument("-v","--verbose",help="Be verbose, say everything program does. Default is False",required=False,action="store_true")
    parser.add_argument("--filename",type=str,help="Name of the measurement set for which weights want to be calculated",required=True,nargs="+")
    parser.add_argument("--ntsol",type=int,help="Solution interval, in timesteps, for your calibration",required=True)
    parser.add_argument("--nfreqsol",type=int,help="Frequency interval, in channels, for your calibration. Default is 8.",required=False,default=8)
    parser.add_argument("--colname",type=str,help="Name of the weights column name you want to save the weights to. Default is CAL_WEIGHT.",required=False,default="CAL_WEIGHT")
    parser.add_argument("--gainfile",type=str,help="Name of the gain file you want to read to rebuild the calibration quality weights."+\
                        " If no file is given, equivalent to rebuilding weights for phase-only calibration.",required=False,default="")
    parser.add_argument("--uvcutkm",type=float,nargs=2,default=[0,2000],required=False,help="uvcut used during calibration, in km.")
    parser.add_argument("--phaseonly",help="Use if calibration was phase-only; this means that gain information doesn't need to be read.",required=False,action="store_true")
    parser.add_argument("--normalise",help="Normalise gains to avoid suppressing long baselines",required=False,action="store_true")
    args=parser.parse_args()
    return vars(args)



### if program is called as main ###
if __name__=="__main__":
    start_time=time.time()
    args        = readArguments()
    mslist      = args["filename"]
    ntsol       = args["ntsol"]
    nfreqsol    = args["nfreqsol"]
    colname     = args["colname"]
    gainfile    = args["gainfile"]
    uvcut       = args["uvcutkm"]
    phaseonly   = args["phaseonly"]
    normalise   = args["normalise"]
    for msname in mslist:
        print "Finding time-covariance weights for: %s"%msname
        covweights=CovWeights(MSName=msname,ntsol=ntsol,gainfile=gainfile,uvcut=uvcut,phaseonly=phaseonly,norm=normalise)
        coefficients=covweights.FindWeights(tcorr=0,colname=colname)
        covweights.SaveWeights(coefficients,colname=colname,AverageOverChannels=True,tcorr=0)
        print "Total runtime: %f min"%((time.time()-start_time)/60.)
