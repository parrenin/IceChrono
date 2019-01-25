#TODO: extend the chronology down to the bedrock by extrapolating the accumulation
#TODO: optinally use a restart file to have a bootstrap method
#TODO: is there an elegant way to unpack the variables vector in the model function?
#TODO: allow to save the correction vector to be able to restart while changing the resolution
#TODO: include some checks for when dDdepth/dz>1
#TODO: Delta-depth observations should be lognormal?
#TODO: we should superpose two charts for ice and air ages, one for the age and one for the uncertainty, since the min age is not always near 0.
#TODO: also compute the prior uncertainties and show them in the figures.
#TODO: the reading of observations does not work if there is only one observation (since the readed matrix is 1D in this case).
#TODO: is there really a computation gain with the change of variable for the correction functions? Avoiding this change of variables would make the code easier to understand. I think there is no gain since solving A^-1 b when we have the LU factorisation of A does not cost more than computing A^-1 * b when we have computed A^-1.


def interp_lin_aver(xp, x, y):
    yp=np.nan*np.zeros(np.size(xp)-1)
    if xp[0]<min(x):
        xmod=np.concatenate((np.array([xp[0]]),x))
        ymod=np.concatenate((np.array([y[0]]),y))
    else:
        xmod=x+0
        ymod=y+0
    if xp[-1]>max(x):
        xmod=np.concatenate((xmod,np.array([xp[-1]])))
        ymod=np.concatenate((ymod,np.array([y[-1]])))
    for i in range(np.size(xp)-1):
        xx=xmod[np.where(np.logical_and(xmod>xp[i],xmod<xp[i+1]))]
        xx=np.concatenate((np.array([xp[i]]),xx,np.array([xp[i+1]])))
        yy=np.interp(xx, xmod, ymod)
        yp[i]=np.sum((yy[1:]+yy[:-1])/2*(xx[1:]-xx[:-1]))/(xp[i+1]-xp[i])
    return yp

def interp_stair_aver(xp, x, y):
    xmod=x+0
    ymod=y+0
    if xp[0]<x[0]:
        xmod=np.concatenate((np.array([xp[0]]),xmod))
        ymod=np.concatenate((np.array([y[0]]),ymod))
    if xp[-1]>x[-1]:
        xmod=np.concatenate((xmod,np.array([xp[-1]])))
        ymod=np.concatenate((ymod,np.array([y[-1]])))
    yint=np.cumsum(np.concatenate((np.array([0]),ymod[:-1]*(xmod[1:]-xmod[:-1]))))
    yp=(np.interp(xp[1:], xmod, yint)-np.interp(xp[:-1], xmod, yint))/(xp[1:]-xp[:-1])     #Maybe this is suboptimal since we compute twice g(xp[i])
    return yp


def gaussian(x):
    return np.exp(-x**2/2)

class Drilling:

    def __init__(self, dlabel):
        self.label=dlabel

    def init(self):

#        print 'Initialization of drilling '+self.label

	self.archive='icecore'
        self.accu_prior_rep='staircase'

        execfile(datadir+'/parameters-AllDrillings.py')
        execfile(datadir+self.label+'/parameters.py')

        self.depth_mid=(self.depth[1:]+self.depth[:-1])/2
        self.depth_inter=(self.depth[1:]-self.depth[:-1])

## We set up the raw model

        if self.calc_a:
            readarray=np.loadtxt(datadir+self.label+'/isotopes.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.iso_depth=readarray[:,0]
            if self.calc_a_method=='fullcorr':
                self.iso_d18Oice=readarray[:,1]
                self.d18Oice=interp.stair_aver(self.depth, self.iso_depth, self.iso_d18Oice)
                self.iso_deutice=readarray[:,2]
                self.deutice=interp_stair_aver(self.depth, self.iso_depth, self.iso_deutice)
                self.iso_d18Osw=readarray[:,3]
                self.d18Osw=interp.stair_aver(self.depth, self.iso_depth, self.iso_d18Osw)
                self.excess=self.deutice-8*self.d18Oice   # dans Uemura : d=excess
                self.a=np.empty_like(self.deutice)
                self.d18Oice_corr=self.d18Oice-self.d18Osw*(1+self.d18Oice/1000)/(1+self.d18Osw/1000)	#Uemura (1)
                self.deutice_corr=self.deutice-8*self.d18Osw*(1+self.deutice/1000)/(1+8*self.d18Osw/1000)	#Uemura et al. (CP, 2012) (2) 
                self.excess_corr=self.deutice_corr-8*self.d18Oice_corr
                self.deutice_fullcorr=self.deutice_corr+self.gamma_source/self.beta_source*self.excess_corr
            elif self.calc_a_method=='deut':
                self.iso_deutice=readarray[:,1]
                self.deutice_fullcorr=interp_stair_aver(self.depth, self.iso_depth, self.iso_deutice)
            elif selc.calc_a_method=='d18O':
                self.d18Oice=readarray[:,1]
                self.deutice_fullcorr=8*interp_stair_aver(self.depth, self.iso_depth, self.iso_d18Oice)
            else:
                print 'Accumulation method not recognized'
                quit()
        else:
            readarray=np.loadtxt(datadir+self.label+'/accu-prior.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.a_depth=readarray[:,0]
            self.a_a=readarray[:,1]
            if readarray.shape[1]>=3:
                self.a_sigma=readarray[:,2]
            if self.accu_prior_rep=='staircase':
                self.a_model=interp_stair_aver(self.depth, self.a_depth, self.a_a)
            elif self.accu_prior_rep=='linear':
                self.a_model=interp_lin_aver(self.depth, self.a_depth, self.a_a)
            else:
                print 'Representation of prior accu scenario not recognized'
            self.a=self.a_model


        
        self.age=np.empty_like(self.depth)
        self.airage=np.empty_like(self.depth)
        
        if self.archive=='icecore':

            readarray=np.loadtxt(datadir+self.label+'/density-prior.txt')
            #        self.density_depth=readarray[:,0]
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.D_depth=readarray[:,0]
            self.D_D=readarray[:,1]
            self.D=np.interp(self.depth_mid, self.D_depth, self.D_D)
            self.iedepth=np.cumsum(np.concatenate((np.array([0]), self.D*self.depth_inter)))
            self.iedepth_mid=(self.iedepth[1:]+self.iedepth[:-1])/2

            if self.calc_tau:
                self.thickness_ie=self.thickness-self.depth[-1]+self.iedepth[-1]
            
            if self.calc_LID:
                if self.depth[0]<self.LID_value:
                    self.LID_depth=np.array([self.depth[0], self.LID_value, self.depth[-1]])
                    self.LID_LID=np.array([self.depth[0], self.LID_value, self.LID_value])
                else:
                    self.LID_depth=np.array([self.depth[0], self.depth[-1]])
                    self.LID_LID=np.array([self.LID_value, self.LID_value])
            else:
    #            self.LID_model=np.loadtxt(datadir+self.label+'/LID-prior.txt')
                readarray=np.loadtxt(datadir+self.label+'/LID-prior.txt')
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.LID_depth=readarray[:,0]
                self.LID_LID=readarray[:,1]
                if readarray.shape[1]>=3:
                    self.LID_sigma=readarray[:,2]
            self.LID_model=np.interp(self.depth, self.LID_depth, self.LID_LID)

            self.Ddepth=np.empty_like(self.depth)
            self.udepth=np.empty_like(self.depth)

#        print 'depth_mid ', np.size(self.depth_mid)
#        print 'zeta ', np.size(self.zeta)
            if self.calc_tau:
                self.thicknessie=self.thickness-self.depth[-1]+self.iedepth[-1]
                self.zeta=(self.thicknessie-self.iedepth_mid)/self.thicknessie  #FIXME: maybe we should use iedepth and thickness_ie here?
                self.tau=np.empty_like(self.depth_mid)
            else:
                readarray=np.loadtxt(datadir+self.label+'/thinning-prior.txt')
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.tau_depth=readarray[:,0]
                self.tau_tau=readarray[:,1]
                if readarray.shape[1]>=3:
                    self.tau_sigma=readarray[:,2]
                self.tau_model=np.interp(self.depth_mid, self.tau_depth, self.tau_tau)
                self.tau=self.tau_model

        self.raw_model()

## Now we set up the correction functions

        if self.start=='restart':
            self.variables=np.loadtxt(datadir+self.label+'/restart.txt')
        elif self.start=='default':
            self.corr_a=np.zeros(np.size(self.corr_a_age))
            if self.archive=='icecore':
                self.corr_LID=np.zeros(np.size(self.corr_LID_age))
                self.corr_tau=np.zeros(np.size(self.corr_tau_depth))
        elif self.start=='random':
            self.corr_a=np.random.normal(loc=0., scale=1., size=np.size(self.corr_a_age))
            if self.archive=='icecore':
                self.corr_LID=np.random.normal(loc=0., scale=1., size=np.size(self.corr_LID_age))
                self.corr_tau=np.random.normal(loc=0., scale=1., size=np.size(self.corr_tau_depth))
        else:
            print 'Start option not recognized.'

## Now we set up the correlation matrices

        self.correlation_corr_a=np.diag(np.ones(np.size(self.corr_a)))
        self.chol_a=np.diag(np.ones(np.size(self.corr_a)))
        if self.archive=='icecore':
            self.correlation_corr_LID=np.diag(np.ones(np.size(self.corr_LID)))
            self.correlation_corr_tau=np.diag(np.ones(np.size(self.corr_tau)))
            self.chol_LID=np.diag(np.ones(np.size(self.corr_LID)))
            self.chol_tau=np.diag(np.ones(np.size(self.corr_tau)))



## Definition of the covariance matrix of the background

        try:
            self.sigmap_corr_a=np.interp(self.corr_a_age, self.fct_age_model(self.a_depth), self.a_sigma)           #FIXME: we should average here since it would be more representative
        except AttributeError:
            print 'Sigma on prior accu scenario not defined in the accu-prior.txt file'

        if self.archive=='icecore':
            try:
                self.sigmap_corr_LID=np.interp(self.corr_LID_age, self.fct_airage_model(self.LID_depth) , self.LID_sigma)          #FIXME: we should average here since it would be more representative
            except AttributeError:
                print 'Sigma on prior LID scenario not defined in the LID-prior.txt file'

            try:
                self.sigmap_corr_tau=np.interp(self.corr_tau_depth, self.tau_depth, self.tau_sigma)           #FIXME: we should average here since it would be more representative
            except AttributeError:
                print 'Sigma on prior thinning scenario not defined in the thinning-prior.txt file'

        self.correlation_corr_a_before=self.correlation_corr_a+0
        if self.archive=='icecore':
            self.correlation_corr_LID_before=self.correlation_corr_LID+0
            self.correlation_corr_tau_before=self.correlation_corr_tau+0

    	filename=datadir+self.label+'/parameters-CovariancePrior-init.py'
        if os.path.isfile(filename):
            execfile(filename)
        else:
            filename=datadir+'/parameters-CovariancePrior-AllDrillings-init.py'
            if os.path.isfile(filename):
                execfile(filename)


        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=cholesky(self.correlation_corr_tau)


        self.variables=np.array([])
#        if self.calc_a==True:
#            self.variables=np.concatenate((self.variables, np.array([self.A0]), np.array([self.beta])))
#        if self.calc_tau==True:
#            self.variables=np.concatenate((self.variables, np.array([self.pprime]), np.array([self.muprime])))
        self.variables=np.concatenate((self.variables, self.corr_tau, self.corr_a, self.corr_LID))


#Reading of observations

        filename=datadir+self.label+'/ice_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.icemarkers_depth=readarray[:,0]
                self.icemarkers_age=readarray[:,1]
                self.icemarkers_sigma=readarray[:,2]
            else:
                self.icemarkers_depth=np.array([])
                self.icemarkers_age=np.array([])
                self.icemarkers_sigma=np.array([])

        filename=datadir+self.label+'/air_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.airmarkers_depth=readarray[:,0]
                self.airmarkers_age=readarray[:,1]
                self.airmarkers_sigma=readarray[:,2]
            else:
                self.airmarkers_depth=np.array([])
                self.airmarkers_age=np.array([])
                self.airmarkers_sigma=np.array([])

        filename=datadir+self.label+'/ice_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.iceintervals_depthtop=readarray[:,0]
                self.iceintervals_depthbot=readarray[:,1]
                self.iceintervals_duration=readarray[:,2]
                self.iceintervals_sigma=readarray[:,3]
            else:
                self.iceintervals_depthtop=np.array([])
                self.iceintervals_depthbot=np.array([])
                self.iceintervals_duration=np.array([])
                self.iceintervals_sigma=np.array([])

        filename=datadir+self.label+'/air_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.airintervals_depthtop=readarray[:,0]
                self.airintervals_depthbot=readarray[:,1]
                self.airintervals_duration=readarray[:,2]
                self.airintervals_sigma=readarray[:,3]
            else:
                self.airintervals_depthtop=np.array([])
                self.airintervals_depthbot=np.array([])
                self.airintervals_duration=np.array([])
                self.airintervals_sigma=np.array([])

        filename=datadir+self.label+'/Ddepth.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.Ddepth_depth=readarray[:,0]
                self.Ddepth_Ddepth=readarray[:,1]
                self.Ddepth_sigma=readarray[:,2]
            else:
                self.Ddepth_depth=np.array([])
                self.Ddepth_Ddepth=np.array([])
                self.Ddepth_sigma=np.array([])


        self.icemarkers_correlation=np.diag(np.ones(np.size(self.icemarkers_depth)))
        self.airmarkers_correlation=np.diag(np.ones(np.size(self.airmarkers_depth)))
        self.iceintervals_correlation=np.diag(np.ones(np.size(self.iceintervals_depthtop)))
        self.airintervals_correlation=np.diag(np.ones(np.size(self.airintervals_depthtop)))
        self.Ddepth_correlation=np.diag(np.ones(np.size(self.Ddepth_depth)))
#        print self.icemarkers_correlation
        filename=datadir+'/parameters-CovarianceObservations-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.icemarkers_depth)>0:
            self.icemarkers_chol=cholesky(self.icemarkers_correlation)
            self.icemarkers_lu_piv=scipy.linalg.lu_factor(np.transpose(self.icemarkers_chol))  #FIXME: we LU factor a triangular matrix. This is suboptimal. We should set lu_piv directly instead.
        if np.size(self.airmarkers_depth)>0:
            self.airmarkers_chol=cholesky(self.airmarkers_correlation)
            self.airmarkers_lu_piv=scipy.linalg.lu_factor(np.transpose(self.airmarkers_chol))
        if np.size(self.iceintervals_depthtop)>0:
            self.iceintervals_chol=cholesky(self.iceintervals_correlation)
            self.iceintervals_lu_piv=scipy.linalg.lu_factor(np.transpose(self.iceintervals_chol))
        if np.size(self.airintervals_depthtop)>0:
            self.airintervals_chol=cholesky(self.airintervals_correlation)
            self.airintervals_lu_piv=scipy.linalg.lu_factor(np.transpose(self.airintervals_chol))
        if np.size(self.Ddepth_depth)>0:
            self.Ddepth_chol=cholesky(self.Ddepth_correlation)
            self.Ddepth_lu_piv=scipy.linalg.lu_factor(np.transpose(self.Ddepth_chol))


    def raw_model(self):



        #Accumulation
        if self.calc_a:
            self.a_model=self.A0*np.exp(self.beta*(self.deutice_fullcorr-self.deutice_fullcorr[0])) #Parrenin et al. (CP, 2007a) 2.3 (6)

        #Thinning
        if self.calc_tau:
            self.p=-1+m.exp(self.pprime)
            self.mu=m.exp(self.muprime)
#            self.s=m.tanh(self.sprime)
            omega_D=1-(self.p+2)/(self.p+1)*(1-self.zeta)+1/(self.p+1)*(1-self.zeta)**(self.p+2)	#Parrenin et al. (CP, 2007a) 2.2 (3)
            omega=self.s*self.zeta+(1-self.s)*omega_D   #Parrenin et al. (CP, 2007a) 2.2 (2)
            self.tau_model=(1-self.mu)*omega+self.mu 

        #udepth
        self.udepth_model=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model*self.depth_inter)))
        
        self.LIDIE_model=self.LID_model*self.Dfirn
        self.ULIDIE_model=np.interp(self.LIDIE_model, self.iedepth, self.udepth_model)

        #Ice age
        self.icelayerthick_model=self.tau_model*self.a_model/self.D
        self.age_model=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model/self.a_model*self.depth_inter)))
            

        #air age
#        self.ice_equiv_depth_model=i_model(np.where(self.udepth_model-self.ULIDIE_model>self.udepth_top, self.udepth_model-self.ULIDIE_model, np.nan))  
        self.ice_equiv_depth_model=np.interp(self.udepth_model-self.ULIDIE_model, self.udepth_model, self.depth)  
        self.Ddepth_model=self.depth-self.ice_equiv_depth_model
        self.airage_model=np.interp(self.ice_equiv_depth_model, self.depth, self.age_model, left=np.nan, right=np.nan)
        self.airlayerthick_model=1/np.diff(self.airage_model)

    def corrected_model(self):

        self.correlation_corr_a_before=self.correlation_corr_a+0
        self.correlation_corr_LID_before=self.correlation_corr_LID+0
        self.correlation_corr_tau_before=self.correlation_corr_tau+0

        filename=datadir+'/parameters-CovariancePrior-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovariancePrior.py'
        if os.path.isfile(filename):
            execfile(filename)

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=cholesky(self.correlation_corr_tau)


        #Accu
        corr=np.dot(self.chol_a,self.corr_a)*self.sigmap_corr_a
        self.a=self.a_model*np.exp(np.interp(self.age_model[:-1], self.corr_a_age, corr)) #FIXME: we should use mid-age and not age

        #Thinning
        self.tau=self.tau_model*np.exp(np.interp(self.depth_mid, self.corr_tau_depth, np.dot(self.chol_tau,self.corr_tau)*self.sigmap_corr_tau))
        self.udepth=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau*self.depth_inter)))
        corr=np.dot(self.chol_LID,self.corr_LID)*self.sigmap_corr_LID
        self.LID=self.LID_model*np.exp(np.interp(self.age_model, self.corr_LID_age, corr))
        self.LIDIE=self.LID*self.Dfirn
        self.ULIDIE=np.interp(self.LIDIE, self.iedepth, self.udepth)

        #Ice age
        self.icelayerthick=self.tau*self.a/self.D
        self.age=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau/self.a*self.depth_inter)))

        self.ice_equiv_depth=np.interp(self.udepth-self.ULIDIE, self.udepth, self.depth)
        self.Ddepth=self.depth-self.ice_equiv_depth
        self.airage=np.interp(self.ice_equiv_depth, self.depth,self.age, left=np.nan, right=np.nan)
        self.airlayerthick=1/np.diff(self.airage)


    def model(self, variables):
        index=0
#        if self.calc_a==True:
#            self.A0=variables[index]
#            self.beta=variables[index+1]
#            index=index+2
#        if self.calc_tau==True:
##            self.p=-1+m.exp(variables[index])
##            self.s=variables[index+1]
##            self.mu=variables[index+2]
##            index=index+3
#            self.pprime=variables[index]
#            self.muprime=variables[index+1]
#            index=index+2
        self.corr_tau=variables[index:index+np.size(self.corr_tau)]
        self.corr_a=variables[index+np.size(self.corr_tau):index+np.size(self.corr_tau)+np.size(self.corr_a)]
        self.corr_LID=variables[index+np.size(self.corr_tau)+np.size(self.corr_a):index+np.size(self.corr_tau)+np.size(self.corr_a)+np.size(self.corr_LID)]

        ##Raw model

        self.raw_model()

        ##Corrected model

        self.corrected_model()

        return np.concatenate((self.age,self.airage,self.Ddepth,self.a,self.tau,self.LID,self.icelayerthick,self.airlayerthick)) 


    def write_init(self):
        self.a_init=self.a
        self.LID_init=self.LID
        self.tau_init=self.tau
        self.icelayerthick_init=self.icelayerthick
        self.airlayerthick_init=self.airlayerthick
        self.age_init=self.age
        self.airage_init=self.airage
        self.Ddepth_init=self.Ddepth

    def fct_age(self, depth):
        return np.interp(depth, self.depth, self.age)

    def fct_age_init(self, depth):
        return np.interp(depth, self.depth, self.age_init)
   
    def fct_age_model(self, depth):
        return np.interp(depth, self.depth,self.age_model)
   
    def fct_airage(self, depth):
        return np.interp(depth, self.depth, self.airage)

    def fct_airage_init(self, depth):
        return np.interp(depth, self.depth, self.airage_init)

    def fct_airage_model(self, depth):
        return np.interp(depth, self.depth, self.airage_model)

    def fct_Ddepth(self, depth):
        return np.interp(depth, self.depth, self.Ddepth)

    def residuals(self, variables):
        self.model(variables)
        resi_corr_a=self.corr_a
        resi_corr_LID=self.corr_LID
        resi_corr_tau=self.corr_tau
        resi_age=(self.fct_age(self.icemarkers_depth)-self.icemarkers_age)/self.icemarkers_sigma
        if np.size(self.icemarkers_depth)>0:
            resi_age=scipy.linalg.lu_solve(self.icemarkers_lu_piv,resi_age)
        resi_airage=(self.fct_airage(self.airmarkers_depth)-self.airmarkers_age)/self.airmarkers_sigma
        if np.size(self.airmarkers_depth)>0:
            resi_airage=scipy.linalg.lu_solve(self.airmarkers_lu_piv,resi_airage)
        resi_iceint=(self.fct_age(self.iceintervals_depthbot)-self.fct_age(self.iceintervals_depthtop)-self.iceintervals_duration)/self.iceintervals_sigma
        if np.size(self.iceintervals_depthtop)>0:
            resi_iceint=scipy.linalg.lu_solve(self.iceintervals_lu_piv,resi_iceint)
        resi_airint=(self.fct_airage(self.airintervals_depthbot)-self.fct_airage(self.airintervals_depthtop)-self.airintervals_duration)/self.airintervals_sigma
        if np.size(self.airintervals_depthtop)>0:
            resi_airint=scipy.linalg.lu_solve(self.airintervals_lu_piv,resi_airint)
        resi_Ddepth=(self.fct_Ddepth(self.Ddepth_depth)-self.Ddepth_Ddepth)/self.Ddepth_sigma
        if np.size(self.Ddepth_depth)>0:
            resi_Ddepth=scipy.linalg.lu_solve(self.Ddepth_lu_piv,resi_Ddepth)
        return np.concatenate((resi_corr_a, resi_corr_LID, resi_corr_tau, resi_age,resi_airage, resi_iceint, resi_airint, resi_Ddepth))


    def cost_function(self):
        cost=np.dot(self.residuals,np.transpose(self.residuals))
        return cost

    def jacobian(self):
        epsilon=np.sqrt(np.diag(self.hess))/100000000.
        model0=self.model(self.variables)
        jacob=np.empty((np.size(model0), np.size(self.variables)))
        for i in np.arange(np.size(self.variables)):
            var=self.variables+0
            var[i]=var[i]+epsilon[i]
            model1=self.model(var)
            jacob[:,i]=(model1-model0)/epsilon[i]
        model0=self.model(self.variables)

        return jacob
    
    
    def optimisation(self) : 
        self.variables,self.hess,self.infodict,mesg,ier=leastsq(self.residuals, self.variables, full_output=1)
        print self.variables
        print self.hess
        return self.variables, self.hess
      
        
    def sigma(self):
        jacob=self.jacobian()

        index=0
        c_model=np.dot(jacob[index:index+np.size(self.age),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.age),:])))
        self.sigma_age=np.sqrt(np.diag(c_model))
        index=index+np.size(self.age)
        c_model=np.dot(jacob[index:index+np.size(self.airage),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.airage),:])))
        self.sigma_airage=np.sqrt(np.diag(c_model))
        index=index+np.size(self.airage)
        c_model=np.dot(jacob[index:index+np.size(self.Ddepth),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.Ddepth),:])))
        self.sigma_Ddepth=np.sqrt(np.diag(c_model))
        index=index+np.size(self.Ddepth)
        c_model=np.dot(jacob[index:index+np.size(self.a),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.a),:])))
        self.sigma_a=np.sqrt(np.diag(c_model))
        index=index+np.size(self.a)
        c_model=np.dot(jacob[index:index+np.size(self.tau),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.tau),:])))
        self.sigma_tau=np.sqrt(np.diag(c_model))
        index=index+np.size(self.tau)
        c_model=np.dot(jacob[index:index+np.size(self.LID),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.LID),:])))
        self.sigma_LID=np.sqrt(np.diag(c_model))
        index=index+np.size(self.LID)
        c_model=np.dot(jacob[index:index+np.size(self.icelayerthick),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.icelayerthick),:])))
        self.sigma_icelayerthick=np.sqrt(np.diag(c_model))
        index=index+np.size(self.icelayerthick)
        c_model=np.dot(jacob[index:index+np.size(self.airlayerthick),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.airlayerthick),:])))
        self.sigma_airlayerthick=np.sqrt(np.diag(c_model))


        self.sigma_a_model=np.interp((self.age_model[1:]+self.age_model[:-1])/2, self.corr_a_age, self.sigmap_corr_a)
        self.sigma_LID_model=np.interp(self.age_model, self.corr_LID_age, self.sigmap_corr_LID)
        self.sigma_tau_model=np.interp(self.depth_mid, self.corr_tau_depth, self.sigmap_corr_tau)

    def sigma_zero(self):

        self.sigma_age=np.zeros_like(self.age)
        self.sigma_airage=np.zeros_like(self.airage)
        self.sigma_Ddepth=np.zeros_like(self.Ddepth)
        self.sigma_a=np.zeros_like(self.a)
        self.sigma_tau=np.zeros_like(self.tau)
        self.sigma_LID=np.zeros_like(self.LID)
        self.sigma_icelayerthick=np.zeros_like(self.icelayerthick)
        self.sigma_airlayerthick=np.zeros_like(self.airlayerthick)
        self.sigma_a_model=np.interp((self.age_model[1:]+self.age_model[:-1])/2, self.corr_a_age, self.sigmap_corr_a)
        self.sigma_LID_model=np.interp(self.age_model, self.corr_LID_age, self.sigmap_corr_LID)
        self.sigma_tau_model=np.interp(self.depth_mid, self.corr_tau_depth, self.sigmap_corr_tau)

        
    

        
    def figures(self):

        mpl.figure(self.label+' thinning')
        mpl.title(self.label+' thinning')
        mpl.xlabel('Thinning')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.tau_init, self.depth_mid, color=color_init, label='Initial')
        mpl.plot(self.tau_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.tau, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.tau-self.sigma_tau, self.tau+self.sigma_tau, color=color_ci)
#        mpl.plot(self.tau+self.sigma_tau, self.depth_mid, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.tau-self.sigma_tau, self.depth_mid, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/thinning.pdf')
        pp.savefig(mpl.figure(self.label+' thinning'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice layer thickness')
        mpl.title(self.label+' ice layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.icelayerthick_init, self.depth_mid, color=color_init, label='Initial')
#        for i in range(np.size(self.iceintervals_duration)):
#            y1=self.iceintervals_depthtop[i]
#            y2=self.iceintervals_depthbot[i]
#            x1=(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
#            x2=(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
#            yserie=np.array([y1,y1,y2,y2,y1])
#            xserie=np.array([x1,x2,x2,x1,x1])
#            if i==0:
#                mpl.plot(xserie,yserie, color=color_obs, label="observations")
#            else:
#                mpl.plot(xserie,yserie, color=color_obs)
        mpl.plot(self.icelayerthick_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.icelayerthick, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.icelayerthick-self.sigma_icelayerthick, self.icelayerthick+self.sigma_icelayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((0,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/icelayerthick.pdf')
        pp.savefig(mpl.figure(self.label+' ice layer thickness'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air layer thickness')
        mpl.title(self.label+' air layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.airlayerthick_init, self.depth_mid, color=color_init, label='Initial')
#        for i in range(np.size(self.airintervals_duration)):
#            y1=self.airintervals_depthtop[i]
#            y2=self.airintervals_depthbot[i]
#            x1=(y2-y1)/(self.airintervals_duration[i]+self.airintervals_sigma[i])
#            x2=(y2-y1)/(self.airintervals_duration[i]-self.airintervals_sigma[i])
#            yserie=np.array([y1,y1,y2,y2,y1])
#            xserie=np.array([x1,x2,x2,x1,x1])
#            if i==0:
#                mpl.plot(xserie,yserie, color=color_obs, label='observations')
#            else:
#                mpl.plot(xserie,yserie, color=color_obs)
        mpl.plot(self.airlayerthick_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.airlayerthick, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.airlayerthick-self.sigma_airlayerthick, self.airlayerthick+self.sigma_airlayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((0, 2*max(self.icelayerthick),self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/airlayerthick.pdf')
        if show_airlayerthick:
            pp.savefig(mpl.figure(self.label+' air layer thickness'))  #Fixme: buggy line on anaconda
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' accumulation')
        mpl.title(self.label+' accumulation')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('Accumulation (m/yr)')
        if show_initial:
            mpl.step(self.age, np.concatenate((self.a_init, np.array([self.a_init[-1]]))), color=color_init, where='post', label='Initial')
        mpl.step(self.age, np.concatenate((self.a_model, np.array([self.a_model[-1]]))), color=color_mod, where='post', label='Prior')
        mpl.step(self.age, np.concatenate((self.a, np.array([self.a[-1]]))), color=color_opt, where='post', label='Posterior +/-$\sigma$')
        mpl.fill_between(self.age[:-1], self.a-self.sigma_a, self.a+self.sigma_a, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/accumulation.pdf')
        pp.savefig(mpl.figure(self.label+' accumulation'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' LID')
        mpl.title(self.label+' LID')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('LID')
        if show_initial:
            mpl.plot(self.age, self.LID_init, color=color_init, label='Initial')
        mpl.plot(self.age, self.LID_model, color=color_mod, label='Prior')
        mpl.plot(self.age, self.LID, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_between(self.age, self.LID-self.sigma_LID, self.LID+self.sigma_LID, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/LID.pdf')
        pp.savefig(mpl.figure(self.label+' LID'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice age')
        mpl.title(self.label+' ice age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.age_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.icemarkers_depth)>0):
            mpl.errorbar(self.icemarkers_age, self.icemarkers_depth, color=color_obs, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="dated horizons")
#        mpl.ylim(mpl.ylim()[::-1])
        for i in range(np.size(self.iceintervals_duration)):
            y1=self.iceintervals_depthtop[i]
            y2=self.iceintervals_depthbot[i]
            x1=self.fct_age(y1)  #(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
            x2=x1+self.iceintervals_duration[i]  #(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
            xseries=np.array([x1,x2,x2,x1,x1])
            yseries=np.array([y1,y1,y2,y2,y1])
            if i==0:
                mpl.plot(xseries, yseries, color=color_di, label="dated intervals")
                mpl.errorbar(x2, y2, color=color_di, xerr=self.iceintervals_sigma[i], capsize=1)
            else:
                mpl.plot(xseries, yseries, color=color_di)
                mpl.errorbar(x2, y2, color=color_di, xerr=self.iceintervals_sigma[i], capsize=1)
#            mpl.arrow(x1, y1, x2-x1, y2-y1, fc=color_di, ec=color_di, head_width=0.02, head_length=0.05)        
#        if (np.size(self.iceintervals_depthtop)>0):
#            mpl.errorbar(self.fct_age(self.iceintervals_depthtop)+self.iceintervals_duration, self.iceintervals_depthbot, color=color_di, xerr=self.iceintervals_sigma, linestyle='', marker='o', markersize='2', label="dated intervals")
        mpl.plot(self.age_model, self.depth, color=color_mod, label='Prior')
        mpl.plot(self.age, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.age-self.sigma_age, self.age+self.sigma_age , color=color_ci)
#        mpl.plot(self.age-self.sigma_age, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_age*scale_ageci, self.depth, color=color_sigma, label='$\sigma$ x'+str(scale_ageci))   
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,self.depth[-1],self.depth[0]))    
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice_age.pdf')
        pp.savefig(mpl.figure(self.label+' ice age'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air age')
        mpl.title(self.label+' air age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.airage_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.airmarkers_depth)>0):
            mpl.errorbar(self.airmarkers_age, self.airmarkers_depth, color=color_obs, xerr=self.airmarkers_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        for i in range(np.size(self.airintervals_duration)):
            y1=self.airintervals_depthtop[i]
            y2=self.airintervals_depthbot[i]
            x1=self.fct_airage(y1)  #(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
            x2=x1+self.airintervals_duration[i]  #(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
            xseries=np.array([x1,x2,x2,x1,x1])
            yseries=np.array([y1,y1,y2,y2,y1])
            if i==0:
                mpl.plot(xseries, yseries, color=color_di, label="dated intervals")
                mpl.errorbar(x2, y2, color=color_di, xerr=self.airintervals_sigma[i], capsize=1)
            else:
                mpl.plot(xseries, yseries, color=color_di)
                mpl.errorbar(x2, y2, color=color_di, xerr=self.airintervals_sigma[i], capsize=1)
        mpl.plot(self.airage_model, self.depth, color=color_mod, label='Prior')
        mpl.fill_betweenx(self.depth, self.airage-self.sigma_airage, self.airage+self.sigma_airage , color=color_ci)
        mpl.plot(self.airage, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
#        mpl.plot(self.airage+self.sigma_airage, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.airage-self.sigma_airage, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_airage*scale_ageci, self.depth, color=color_sigma, label='$\sigma$ x'+str(scale_ageci))  
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,self.depth[-1],self.depth[0]))    
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air_age.pdf')
        pp.savefig(mpl.figure(self.label+' air age'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' Ddepth')
        mpl.title(self.label+' $\Delta$depth')
        mpl.xlabel('$\Delta$depth (m)')
        mpl.ylabel('Air depth (m)')
        if show_initial:
            mpl.plot(self.Ddepth_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.Ddepth_depth)>0):
            mpl.errorbar(self.Ddepth_Ddepth, self.Ddepth_depth, color=color_obs, xerr=self.Ddepth_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        mpl.plot(self.Ddepth_model, self.depth, color=color_mod, label='Prior')
        mpl.plot(self.Ddepth, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.Ddepth-self.sigma_Ddepth, self.Ddepth+self.sigma_Ddepth, color=color_ci)
#        mpl.plot(self.Ddepth+self.sigma_Ddepth, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.Ddepth-self.sigma_Ddepth, self.depth, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/Ddepth.pdf')
        pp.savefig(mpl.figure(self.label+' Ddepth'))
        pp.close()
        if not show_figures:
            mpl.close()


    def save(self):
        output=np.vstack((self.depth,self.age,self.sigma_age,self.airage,self.sigma_airage,np.append(self.a,self.a[-1]),np.append(self.sigma_a,self.sigma_a[-1]),np.append(self.tau,self.tau[-1]),np.append(self.sigma_tau,self.sigma_tau[-1]),self.LID,self.sigma_LID, self.Ddepth,self.sigma_Ddepth,np.append(self.a_model,self.a_model[-1]),np.append(self.sigma_a_model,self.sigma_a_model[-1]),np.append(self.tau_model,self.tau_model[-1]),np.append(self.sigma_tau_model,self.sigma_tau_model[-1]),self.LID_model,self.sigma_LID_model,np.append(self.icelayerthick,self.icelayerthick[-1]),np.append(self.sigma_icelayerthick,self.sigma_icelayerthick[-1]),np.append(self.airlayerthick,self.airlayerthick[-1]),np.append(self.sigma_airlayerthick,self.sigma_airlayerthick[-1])))
        with open(datadir+self.label+'/output.txt','w') as f:
            f.write('#depth\tage\tsigma_age\tair_age\tsigma_air_age\taccu\tsigma_accu\tthinning\tsigma_thinning\tLID\tsigma_LID\tDdepth\tsigma_Ddepth\taccu_model\tsigma_accu_model\tthinning_model\tsigma_thinning_model\tLID_model\tsigma_LID_model\ticelayerthick\tsigma_icelayerthick\tairlayerthick\tsigma_airlayerthick\n')
            np.savetxt(f,np.transpose(output), delimiter='\t')
        np.savetxt(datadir+self.label+'/restart.txt',np.transpose(self.variables))
    
#    def udepth_save(self):
#        np.savetxt(datadir+self.label+'/udepth.txt',self.udepth)


class DrillingPair:

    def __init__(self, D1, D2):
        self.D1=D1
        self.D2=D2


    def init(self):
        self.label=self.D1.label+'-'+self.D2.label
#        print 'Initialization of drilling pair ',self.label


#TODO: allow to have either dlabel1+'-'dlabel2 or dlbel2+'-'dlabel1 as directory
        filename=datadir+self.D1.label+'-'+self.D2.label+'/ice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.iceicemarkers_depth1=readarray[:,0]
            self.iceicemarkers_depth2=readarray[:,1]
            self.iceicemarkers_sigma=readarray[:,2]
        else:
            self.iceicemarkers_depth1=np.array([])
            self.iceicemarkers_depth2=np.array([])
            self.iceicemarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/air_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.airairmarkers_depth1=readarray[:,0]
            self.airairmarkers_depth2=readarray[:,1]
            self.airairmarkers_sigma=readarray[:,2]
        else:
            self.airairmarkers_depth1=np.array([])
            self.airairmarkers_depth2=np.array([])
            self.airairmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/iceair_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.iceairmarkers_depth1=readarray[:,0]
            self.iceairmarkers_depth2=readarray[:,1]
            self.iceairmarkers_sigma=readarray[:,2]
        else:
            self.iceairmarkers_depth1=np.array([])
            self.iceairmarkers_depth2=np.array([])
            self.iceairmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/airice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.airicemarkers_depth1=readarray[:,0]
            self.airicemarkers_depth2=readarray[:,1]
            self.airicemarkers_sigma=readarray[:,2]
        else:
            self.airicemarkers_depth1=np.array([])
            self.airicemarkers_depth2=np.array([])
            self.airicemarkers_sigma=np.array([])


        self.iceicemarkers_correlation=np.diag(np.ones(np.size(self.iceicemarkers_depth1)))
        self.airairmarkers_correlation=np.diag(np.ones(np.size(self.airairmarkers_depth1)))
        self.iceairmarkers_correlation=np.diag(np.ones(np.size(self.iceairmarkers_depth1)))
        self.airicemarkers_correlation=np.diag(np.ones(np.size(self.airicemarkers_depth1)))
        filename=datadir+'/parameters-CovarianceObservations-AllDrillingPairs.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.iceicemarkers_depth1)>0:
            self.iceicemarkers_chol=cholesky(self.iceicemarkers_correlation)
            self.iceicemarkers_lu_piv=scipy.linalg.lu_factor(self.iceicemarkers_chol)
        if np.size(self.airairmarkers_depth1)>0:
            self.airairmarkers_chol=cholesky(self.airairmarkers_correlation)
            self.airairmarkers_lu_piv=scipy.linalg.lu_factor(self.airairmarkers_chol)
        if np.size(self.iceairmarkers_depth1)>0:
            self.iceairmarkers_chol=cholesky(self.iceairmarkers_correlation)
            self.iceairmarkers_lu_piv=scipy.linalg.lu_factor(self.iceairmarkers_chol)
        if np.size(self.airicemarkers_depth1)>0:
            self.airicemarkers_chol=cholesky(self.airicemarkers_correlation)
            self.airicemarkers_lu_piv=scipy.linalg.lu_factor(self.airicemarkers_chol)


    def residuals(self):

        resi_iceice=(self.D1.fct_age(self.iceicemarkers_depth1)-self.D2.fct_age(self.iceicemarkers_depth2))/self.iceicemarkers_sigma
        if np.size(self.iceicemarkers_depth1)>0:
            resi_iceice=scipy.linalg.lu_solve(self.iceicemarkers_lu_piv,resi_iceice)
        resi_airair=(self.D1.fct_airage(self.airairmarkers_depth1)-self.D2.fct_airage(self.airairmarkers_depth2))/self.airairmarkers_sigma
        if np.size(self.airairmarkers_depth1)>0:
            resi_airair=scipy.linalg.lu_solve(self.airairmarkers_lu_piv,resi_airair)
        resi_iceair=(self.D1.fct_age(self.iceairmarkers_depth1)-self.D2.fct_airage(self.iceairmarkers_depth2))/self.iceairmarkers_sigma
        if np.size(self.iceairmarkers_depth1)>0:
            resi_iceair=scipy.linalg.lu_solve(self.iceairmarkers_lu_piv,resi_iceair)
        resi_airice=(self.D1.fct_airage(self.airicemarkers_depth1)-self.D2.fct_age(self.airicemarkers_depth2))/self.airicemarkers_sigma
        if np.size(self.airicemarkers_depth1)>0:
            resi_airice=scipy.linalg.lu_solve(self.airicemarkers_lu_piv,resi_airice)
        resi=np.concatenate((resi_iceice,resi_airair,resi_iceair,resi_airice))
        
        return resi
    

    def figures(self):

        if not os.path.isdir(datadir+self.label):
            os.mkdir(datadir+self.label)


        mpl.figure(self.label+' ice-ice')
        mpl.xlabel(self.D1.label+' ice age (yr b1950)')
        mpl.ylabel(self.D2.label+' ice age (yr b1950)')
        if (np.size(self.iceicemarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_age_init(self.iceicemarkers_depth1),self.D2.fct_age_init(self.iceicemarkers_depth2), color=color_init, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_age_model(self.iceicemarkers_depth1),self.D2.fct_age_model(self.iceicemarkers_depth2), color=color_mod, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_age(self.iceicemarkers_depth1),self.D2.fct_age(self.iceicemarkers_depth2), color=color_opt, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice-ice.pdf')
        pp.savefig(mpl.figure(self.label+' ice-ice'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air-air')
        mpl.xlabel(self.D1.label+' air age (yr b1950)')
        mpl.ylabel(self.D2.label+' air age (yr b1950)')
        if (np.size(self.airairmarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_airage_init(self.airairmarkers_depth1),self.D2.fct_airage_init(self.airairmarkers_depth2), color=color_init, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_airage_model(self.airairmarkers_depth1),self.D2.fct_airage_model(self.airairmarkers_depth2), color=color_mod, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_airage(self.airairmarkers_depth1),self.D2.fct_airage(self.airairmarkers_depth2), color=color_opt, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air-air.pdf')
        pp.savefig(mpl.figure(self.label+' air-air'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice-air')
        mpl.xlabel(self.D1.label+' ice age (yr b1950)')
        mpl.ylabel(self.D2.label+' air age (yr b1950)')
        if (np.size(self.iceairmarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_age_init(self.iceairmarkers_depth1),self.D2.fct_airage_init(self.iceairmarkers_depth2), color=color_init, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_age_model(self.iceairmarkers_depth1),self.D2.fct_airage_model(self.iceairmarkers_depth2), color=color_mod, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_age(self.iceairmarkers_depth1),self.D2.fct_airage(self.iceairmarkers_depth2), color=color_opt, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice-air.pdf')
        pp.savefig(mpl.figure(self.label+' ice-air'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air-ice')
        mpl.xlabel(self.D1.label+' air age (yr b1950)')
        mpl.ylabel(self.D2.label+' ice age (yr b1950)')
        if (np.size(self.airicemarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_airage_init(self.airicemarkers_depth1),self.D2.fct_age_init(self.airicemarkers_depth2), color=color_init, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_airage_model(self.airicemarkers_depth1),self.D2.fct_age_model(self.airicemarkers_depth2), color=color_mod, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_airage(self.airicemarkers_depth1),self.D2.fct_age(self.airicemarkers_depth2), color=color_opt, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air-ice.pdf')
        pp.savefig(mpl.figure(self.label+' air-ice'))
        pp.close()
        if not show_figures:
            mpl.close()

