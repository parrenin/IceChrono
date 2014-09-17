#TODO: extend the chronology down to the bedrock by extrapolating the accumulation
#TODO: optinally use a restart file to have a bootstrap method
#TODO: is there an elegant way to unpack the variables vector in the model function?
#TODO: allow to save the correction vector to be able to restart while changing the resolution
#TODO: include some checks for when dDdepth/dz>1


def interp1d_extrap(x,y):
    def f(xp):
        g=interp1d(x,y, bounds_error=False)
        return np.where(xp<x[0],y[0],np.where(xp>x[-1],y[-1],g(xp)))    
    return f

def interp1d_stair_aver_extrap(x, y):
    def f(xp):
        xmod=x+0
        ymod=y+0
        if xp[0]<x[0]:
            xmod=np.concatenate((np.array([xp[0]]),xmod))
            ymod=np.concatenate((np.array([y[0]]),ymod))
        if xp[-1]>x[-1]:
            xmod=np.concatenate((xmod,np.array([xp[-1]])))
            ymod=np.concatenate((ymod,np.array([y[-1]])))
        yint=np.cumsum(np.concatenate((np.array([0]),ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g=interp1d(xmod,yint)
        yp=(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1])     #Maybe this is suboptimal since we compute twice g(xp[i])
        return yp

    return f


class Drilling:

    def __init__(self, dlabel):
        self.label=dlabel

    def init(self):

#        print 'Initialization of drilling '+self.label

        execfile(datadir+'parameters-AllDrillings.py')
        execfile(datadir+self.label+'/parameters.py')

        self.depth_mid=(self.depth[1:]+self.depth[:-1])/2
        self.depth_inter=(self.depth[1:]-self.depth[:-1])

## We set up the raw model

        if self.calc_a:
            readarray=np.loadtxt(datadir+self.label+'/isotopes.txt')
            self.iso_depth=readarray[:,0]
            if self.calc_a_method=='fullcorr':
                self.iso_d18Oice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Oice)
                self.d18Oice=f(self.depth)
                self.iso_deutice=readarray[:,2]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_deutice)
                self.deutice=f(self.depth)
                self.iso_d18Osw=readarray[:,3]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Osw)
                self.d18Osw=f(self.depth)
                self.excess=self.deutice-8*self.d18Oice   # dans Uemura : d=excess
                self.a=np.empty_like(self.deutice)
                self.d18Oice_corr=self.d18Oice-self.d18Osw*(1+self.d18Oice/1000)/(1+self.d18Osw/1000)	#Uemura (1)
                self.deutice_corr=self.deutice-8*self.d18Osw*(1+self.deutice/1000)/(1+8*self.d18Osw/1000)	#Uemura et al. (CP, 2012) (2) 
                self.excess_corr=self.deutice_corr-8*self.d18Oice_corr
                self.deutice_fullcorr=self.deutice_corr+self.gamma_source/self.beta_source*self.excess_corr
            elif self.calc_a_method=='deut':
                self.iso_deutice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_deutice)
                self.deutice_fullcorr=f(self.depth)
            elif selc.calc_a_method=='d18O':
                self.d18Oice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Oice)
                self.deutice_fullcorr=8*f(depth)
            else:
                print 'Accumulation method not recognized'
                quit()
        else:
            readarray=np.loadtxt(datadir+self.label+'/accu-prior.txt')
            self.a_depth=readarray[:,0]
            self.a_a=readarray[:,1]
            f=interp1d_stair_aver_extrap(self.a_depth, self.a_a)
            self.a_model=f(self.depth)
            self.a=self.a_model

        
        self.age=np.empty_like(self.depth)
        self.gage=np.empty_like(self.depth)
        

        readarray=np.loadtxt(datadir+self.label+'/density-prior.txt')
#        self.density_depth=readarray[:,0]
        self.D_depth=readarray[:,0]
        self.D_D=readarray[:,1]
        f=interp1d_extrap(self.D_depth, self.D_D)
        self.D=f(self.depth_mid)

        self.iedepth=np.cumsum(np.concatenate((np.array([0]), self.D*self.depth_inter)))
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
            self.LID_depth=readarray[:,0]
            self.LID_LID=readarray[:,1]
        f=interp1d_extrap(self.LID_depth, self.LID_LID)
        self.LID_model=f(self.depth)




        self.Ddepth=np.empty_like(self.depth)
        self.udepth=np.empty_like(self.depth)

#        print 'depth_mid ', np.size(self.depth_mid)
        self.zeta=(self.thickness-self.depth_mid)/self.thickness
#        print 'zeta ', np.size(self.zeta)
        if self.calc_tau:
            self.tau=np.empty_like(self.depth_mid)
        else:
            readarray=np.loadtxt(datadir+self.label+'/thinning-prior.txt')
            self.tau_depth=readarray[:,0]
            self.tau_tau=readarray[:,1]
            f=interp1d_extrap(self.tau_depth, self.tau_tau)
            self.tau_model=f(self.depth_mid)
            self.tau=self.tau_model

        self.raw_model()

## Now we set up the correction functions

        self.corr_a=np.zeros(np.size(self.corr_a_age))
        self.corr_LID=np.zeros(np.size(self.corr_LID_age))
        self.corr_tau=np.zeros(np.size(self.corr_tau_depth))

        self.correlation_corr_a=np.empty((np.size(self.corr_a),np.size(self.corr_a)))



        self.correlation_corr_LID=np.empty((np.size(self.corr_LID),np.size(self.corr_LID)))


#        print 'depth ', np.size(self.depth)
#        print 'udepth_init ', np.size(self.udepth_init)

        self.correlation_corr_tau=np.empty((np.size(self.corr_tau),np.size(self.corr_tau)))



## Definition of the covariance matrix of the background

        self.correlation_corr_a_before=self.correlation_corr_a+0
        self.correlation_corr_LID_before=self.correlation_corr_LID+0
        self.correlation_corr_tau_before=self.correlation_corr_tau+0

        filename=datadir+'parameters-CovariancePrior-AllDrillings-init.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'parameters-CovariancePrior-init.py'
        if os.path.isfile(filename):
            execfile(filename)

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=np.linalg.cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=np.linalg.cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=np.linalg.cholesky(self.correlation_corr_tau)


        self.variables=np.array([])
        if self.calc_a==True:
            self.variables=np.concatenate((self.variables, np.array([self.A0]), np.array([self.beta])))
        if self.calc_tau==True:
            self.variables=np.concatenate((self.variables, np.array([self.pprime]), np.array([self.muprime])))
        self.variables=np.concatenate((self.variables, self.corr_tau, self.corr_a, self.corr_LID))

        if self.restart:
            self.variables=np.loadtxt(datadir+self.label+'/restart.txt')


#Reading of observations

        filename=datadir+self.label+'/ice_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                self.icemarkers_depth=readarray[:,0]
                self.icemarkers_age=readarray[:,1]
                self.icemarkers_sigma=readarray[:,2]
            else:
                self.icemarkers_depth=np.array([])
                self.icemarkers_age=np.array([])
                self.icemarkers_sigma=np.array([])

        filename=datadir+self.label+'/gas_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                self.gasmarkers_depth=readarray[:,0]
                self.gasmarkers_age=readarray[:,1]
                self.gasmarkers_sigma=readarray[:,2]
            else:
                self.gasmarkers_depth=np.array([])
                self.gasmarkers_age=np.array([])
                self.gasmarkers_sigma=np.array([])

        filename=datadir+self.label+'/ice_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                self.iceintervals_depthtop=readarray[:,0]
                self.iceintervals_depthbot=readarray[:,1]
                self.iceintervals_duration=readarray[:,2]
                self.iceintervals_sigma=readarray[:,3]
            else:
                self.iceintervals_depthtop=np.array([])
                self.iceintervals_depthbot=np.array([])
                self.iceintervals_duration=np.array([])
                self.iceintervals_sigma=np.array([])

        filename=datadir+self.label+'/gas_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                self.gasintervals_depthtop=readarray[:,0]
                self.gasintervals_depthbot=readarray[:,1]
                self.gasintervals_duration=readarray[:,2]
                self.gasintervals_sigma=readarray[:,3]
            else:
                self.gasintervals_depthtop=np.array([])
                self.gasintervals_depthbot=np.array([])
                self.gasintervals_duration=np.array([])
                self.gasintervals_sigma=np.array([])

        filename=datadir+self.label+'/Ddepth.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                self.Ddepth_depth=readarray[:,0]
                self.Ddepth_Ddepth=readarray[:,1]
                self.Ddepth_sigma=readarray[:,2]
            else:
                self.Ddepth_depth=np.array([])
                self.Ddepth_Ddepth=np.array([])
                self.Ddepth_sigma=np.array([])


        self.icemarkers_correlation=np.diag(np.ones(np.size(self.icemarkers_depth)))
        self.gasmarkers_correlation=np.diag(np.ones(np.size(self.gasmarkers_depth)))
        self.iceintervals_correlation=np.diag(np.ones(np.size(self.iceintervals_depthtop)))
        self.gasintervals_correlation=np.diag(np.ones(np.size(self.gasintervals_depthtop)))
        self.Ddepth_correlation=np.diag(np.ones(np.size(self.Ddepth_depth)))
#        print self.icemarkers_correlation
        filename=datadir+'parameters-CovarianceObservations-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.icemarkers_depth)>0:
            self.icemarkers_chol=np.linalg.cholesky(self.icemarkers_correlation)
            self.icemarkers_lu_piv=scipy.linalg.lu_factor(self.icemarkers_chol)  #FIXME: we LU factor a triangular matrix. This is suboptimal. We should set lu_piv directly instead.
        if np.size(self.gasmarkers_depth)>0:
            self.gasmarkers_chol=np.linalg.cholesky(self.gasmarkers_correlation)
            self.gasmarkers_lu_piv=scipy.linalg.lu_factor(self.gasmarkers_chol)
        if np.size(self.iceintervals_depthtop)>0:
            self.iceintervals_chol=np.linalg.cholesky(self.iceintervals_correlation)
            self.iceintervals_lu_piv=scipy.linalg.lu_factor(self.iceintervals_chol)
        if np.size(self.gasintervals_depthtop)>0:
            self.gasintervals_chol=np.linalg.cholesky(self.gasintervals_correlation)
            self.gasintervals_lu_piv=scipy.linalg.lu_factor(self.gasintervals_chol)
        if np.size(self.Ddepth_depth)>0:
            self.Ddepth_chol=np.linalg.cholesky(self.Ddepth_correlation)
            self.Ddepth_lu_piv=scipy.linalg.lu_factor(self.Ddepth_chol)

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
        
        g_model=interp1d(self.iedepth, self.udepth_model)
        self.LIDIE_model=self.LID_model*self.Dfirn
        self.ULIDIE_model=g_model(self.LIDIE_model)
        i_model=interp1d(self.udepth_model, self.depth)

        #Ice age
        self.icelayerthick_model=self.tau_model*self.a_model/self.D
        self.age_model=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model/self.a_model*self.depth_inter)))
            
        f_model=interp1d(self.depth, self.age_model, bounds_error=False, fill_value=np.nan)

        #gas age
        self.ice_equiv_depth_model=i_model(np.where(self.udepth_model-self.ULIDIE_model>self.udepth_top, self.udepth_model-self.ULIDIE_model, np.nan))  
        self.Ddepth_model=self.depth-self.ice_equiv_depth_model
        self.gage_model=f_model(self.ice_equiv_depth_model)
        self.gaslayerthick_model=1/np.diff(self.gage_model)

    def corrected_model(self):

        self.correlation_corr_a_before=self.correlation_corr_a+0
        self.correlation_corr_LID_before=self.correlation_corr_LID+0
        self.correlation_corr_tau_before=self.correlation_corr_tau+0

        filename=datadir+'parameters-CovariancePrior-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'parameters-CovariancePrior.py'
        if os.path.isfile(filename):
            execfile(filename)

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=np.linalg.cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=np.linalg.cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=np.linalg.cholesky(self.correlation_corr_tau)


        #Accu
        corr=np.dot(self.chol_a,self.corr_a)*self.sigmap_corr_a
        j=interp1d_extrap(self.corr_a_age, corr)
        self.a=self.a_model*np.exp(j(self.age_model[:-1])) #FIXME: we should use mid-age and not age

        #Thinning
        h=interp1d(self.corr_tau_depth, np.dot(self.chol_tau,self.corr_tau)*self.sigmap_corr_tau)
        self.tau=self.tau_model*np.exp(h(self.depth_mid))
        self.udepth=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau*self.depth_inter)))
        g=interp1d(self.iedepth, self.udepth)
        corr=np.dot(self.chol_LID,self.corr_LID)*self.sigmap_corr_LID
        j=interp1d_extrap(self.corr_LID_age, corr)
        self.LID=self.LID_model*np.exp(j(self.age_model))
        self.LIDIE=self.LID*self.Dfirn
        self.ULIDIE=g(self.LIDIE)
        i=interp1d(self.udepth, self.depth)

        #Ice age
        self.icelayerthick=self.tau*self.a/self.D
        self.age=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau/self.a*self.depth_inter)))
        f=interp1d(self.depth,self.age, bounds_error=False, fill_value=np.nan)

        self.ice_equiv_depth=i(np.where(self.udepth-self.LIDIE>self.udepth_top, self.udepth-self.LIDIE, np.nan))
        self.Ddepth=self.depth-self.ice_equiv_depth
        self.gage=f(self.ice_equiv_depth)
        self.gaslayerthick=1/np.diff(self.gage)


    def model(self, variables):
        index=0
        if self.calc_a==True:
            self.A0=variables[index]
            self.beta=variables[index+1]
            index=index+2
        if self.calc_tau==True:
#            self.p=-1+m.exp(variables[index])
#            self.s=variables[index+1]
#            self.mu=variables[index+2]
#            index=index+3
            self.pprime=variables[index]
            self.muprime=variables[index+1]
            index=index+2
        self.corr_tau=variables[index:index+np.size(self.corr_tau)]
        self.corr_a=variables[index+np.size(self.corr_tau):index+np.size(self.corr_tau)+np.size(self.corr_a)]
        self.corr_LID=variables[index+np.size(self.corr_tau)+np.size(self.corr_a):index+np.size(self.corr_tau)+np.size(self.corr_a)+np.size(self.corr_LID)]

        ##Raw model

        self.raw_model()

        ##Corrected model

        self.corrected_model()


        return np.concatenate((self.age,self.gage,self.Ddepth,self.a,self.tau,self.LID,self.icelayerthick,self.gaslayerthick)) 

    def fct_age(self, depth):
        f=interp1d(self.depth,self.age)
        return f(depth)
   
    def fct_age_model(self, depth):
        f=interp1d(self.depth,self.age_model)
        return f(depth)
   
    def fct_gage(self, depth):
        f=interp1d(self.depth,self.gage)
        return f(depth)

    def fct_gage_model(self, depth):
        f=interp1d(self.depth,self.gage_model)
        return f(depth)

    def fct_Ddepth(self, depth):
        f=interp1d(self.depth,self.Ddepth)
        return f(depth)

    def residuals(self, variables):
        self.model(variables)
        resi_corr_a=self.corr_a
        resi_corr_LID=self.corr_LID
        resi_corr_tau=self.corr_tau
        resi_age=(self.fct_age(self.icemarkers_depth)-self.icemarkers_age)/self.icemarkers_sigma
        if np.size(self.icemarkers_depth)>0:
            resi_age=scipy.linalg.lu_solve(self.icemarkers_lu_piv,resi_age)
        resi_gage=(self.fct_gage(self.gasmarkers_depth)-self.gasmarkers_age)/self.gasmarkers_sigma
        if np.size(self.gasmarkers_depth)>0:
            resi_gage=scipy.linalg.lu_solve(self.gasmarkers_lu_piv,resi_gage)
        resi_iceint=(self.fct_age(self.iceintervals_depthbot)-self.fct_age(self.iceintervals_depthtop)-self.iceintervals_duration)/self.iceintervals_sigma
        if np.size(self.iceintervals_depthtop)>0:
            resi_iceint=scipy.linalg.lu_solve(self.iceintervals_lu_piv,resi_iceint)
        resi_gasint=(self.fct_gage(self.gasintervals_depthbot)-self.fct_gage(self.gasintervals_depthtop)-self.gasintervals_duration)/self.gasintervals_sigma
        if np.size(self.gasintervals_depthtop)>0:
            resi_gasint=scipy.linalg.lu_solve(self.gasintervals_lu_piv,resi_gasint)
        resi_Ddepth=(self.fct_Ddepth(self.Ddepth_depth)-self.Ddepth_Ddepth)/self.Ddepth_sigma
        if np.size(self.Ddepth_depth)>0:
            resi_Ddepth=scipy.linalg.lu_solve(self.Ddepth_lu_piv,resi_Ddepth)
        return np.concatenate((resi_corr_a, resi_corr_LID, resi_corr_tau, resi_age,resi_gage, resi_iceint, resi_gasint, resi_Ddepth))


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
        c_model=np.dot(jacob[index:index+np.size(self.gage),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.gage),:])))
        self.sigma_gage=np.sqrt(np.diag(c_model))
        index=index+np.size(self.gage)
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
        c_model=np.dot(jacob[index:index+np.size(self.gaslayerthick),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.gaslayerthick),:])))
        self.sigma_gaslayerthick=np.sqrt(np.diag(c_model))

        f=interp1d_extrap(self.corr_a_age, self.sigmap_corr_a)
        self.sigma_a_model=f((self.age_model[1:]+self.age_model[:-1])/2)
        f=interp1d_extrap(self.corr_LID_age, self.sigmap_corr_LID)
        self.sigma_LID_model=f(self.age_model)
        f=interp1d_extrap(self.corr_tau_depth, self.sigmap_corr_tau)
        self.sigma_tau_model=f(self.depth_mid)

        
    

    def display_init(self):
        
    

        mpl.figure(self.label+' thinning')
        mpl.title(self.label+' thinning')
        mpl.xlabel('Thinning')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.tau, self.depth_mid, color=color_init, label='Initial')

        mpl.figure(self.label+' ice layer thickness')
        mpl.title(self.label+' ice layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.icelayerthick, self.depth_mid, color=color_init, label='Initial')
        for i in range(np.size(self.iceintervals_duration)):
            y1=self.iceintervals_depthtop[i]
            y2=self.iceintervals_depthbot[i]
            x1=(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
            x2=(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
            yserie=np.array([y1,y1,y2,y2,y1])
            xserie=np.array([x1,x2,x2,x1,x1])
            if i==0:
                mpl.plot(xserie,yserie, color=color_obs, label="observations")
            else:
                mpl.plot(xserie,yserie, color=color_obs)

        mpl.figure(self.label+' gas layer thickness')
        mpl.title(self.label+' gas layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.gaslayerthick, self.depth_mid, color=color_init, label='Initial')
        for i in range(np.size(self.gasintervals_duration)):
            y1=self.gasintervals_depthtop[i]
            y2=self.gasintervals_depthbot[i]
            x1=(y2-y1)/(self.gasintervals_duration[i]+self.gasintervals_sigma[i])
            x2=(y2-y1)/(self.gasintervals_duration[i]-self.gasintervals_sigma[i])
            yserie=np.array([y1,y1,y2,y2,y1])
            xserie=np.array([x1,x2,x2,x1,x1])
            if i==0:
                mpl.plot(xserie,yserie, color=color_obs, label='observations')
            else:
                mpl.plot(xserie,yserie, color=color_obs)


        mpl.figure(self.label+' accumulation')
        mpl.title(self.label+' accumulation')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('Accumulation (m/an)')

        mpl.figure(self.label+' LID')
        mpl.title(self.label+' LID')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('LID')

        mpl.figure(self.label+' ice age')
        mpl.title(self.label+' ice age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.age, self.depth, color=color_init, label='Initial')
        mpl.errorbar(self.icemarkers_age, self.icemarkers_depth, color=color_obs, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])

        mpl.figure(self.label+' gas age')
        mpl.title(self.label+' gas age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.gage, self.depth, color=color_init, label='Initial')
        mpl.errorbar(self.gasmarkers_age, self.gasmarkers_depth, color=color_obs, xerr=self.gasmarkers_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        
        mpl.figure(self.label+' Ddepth')
        mpl.title(self.label+' $\Delta$depth')
        mpl.xlabel('$\Delta$depth (m)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.Ddepth, self.depth, color=color_init, label='Initial')
        mpl.errorbar(self.Ddepth_Ddepth, self.Ddepth_depth, color=color_obs, xerr=self.Ddepth_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        
    def display_final(self):


        mpl.figure(self.label+' thinning')
        mpl.plot(self.tau_model, self.depth_mid, color=color_mod, label='Model')
        mpl.plot(self.tau, self.depth_mid, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.tau-self.sigma_tau, self.tau+self.sigma_tau, color=color_ci)
#        mpl.plot(self.tau+self.sigma_tau, self.depth_mid, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.tau-self.sigma_tau, self.depth_mid, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[0],self.depth[-1]))
        mpl.legend(loc=4)
        mpl.ylim(mpl.ylim()[::-1])
        pp=PdfPages(datadir+self.label+'/thinning.pdf')
        pp.savefig(mpl.figure(self.label+' thinning'))
        pp.close()

        mpl.figure(self.label+' ice layer thickness')
        mpl.plot(self.icelayerthick_model, self.depth_mid, color=color_mod, label='Model')
        mpl.plot(self.icelayerthick, self.depth_mid, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.icelayerthick-self.sigma_icelayerthick, self.icelayerthick+self.sigma_icelayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[0],self.depth[-1]))
        mpl.legend(loc=4)
        mpl.ylim(mpl.ylim()[::-1])
        pp=PdfPages(datadir+self.label+'/icelayerthick.pdf')
        pp.savefig(mpl.figure(self.label+' ice layer thickness'))
        pp.close()

        mpl.figure(self.label+' gas layer thickness')
        mpl.plot(self.gaslayerthick_model, self.depth_mid, color=color_mod, label='Model')
        mpl.plot(self.gaslayerthick, self.depth_mid, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.gaslayerthick-self.sigma_gaslayerthick, self.gaslayerthick+self.sigma_gaslayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((0, 2*max(self.icelayerthick),self.depth[0],self.depth[-1]))
        mpl.legend(loc=4)
        mpl.ylim(mpl.ylim()[::-1])
        pp=PdfPages(datadir+self.label+'/gaslayerthick.pdf')
        pp.savefig(mpl.figure(self.label+' gas layer thickness'))
        pp.close()

        mpl.figure(self.label+' accumulation')
        if show_initial:
            mpl.step(self.age, np.concatenate((self.a_init, np.array([self.a_init[-1]]))), color=color_init, where='post', label='Initial')
        mpl.step(self.age, np.concatenate((self.a_model, np.array([self.a_model[-1]]))), color=color_mod, where='post', label='Model')
        mpl.step(self.age, np.concatenate((self.a, np.array([self.a[-1]]))), color=color_opt, where='post', label='Corrected +/-$\sigma$')
        mpl.fill_between(self.age[:-1], self.a-self.sigma_a, self.a+self.sigma_a, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend()
        pp=PdfPages(datadir+self.label+'/accumulation.pdf')
        pp.savefig(mpl.figure(self.label+' accumulation'))
        pp.close()

        mpl.figure(self.label+' LID')
        if show_initial:
            mpl.plot(self.age, self.LID_init, color=color_init, label='Initial')
        mpl.plot(self.age, self.LID_model, color=color_mod, label='Model')
        mpl.plot(self.age, self.LID, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_between(self.age, self.LID-self.sigma_LID, self.LID+self.sigma_LID, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend()
        pp=PdfPages(datadir+self.label+'/LID.pdf')
        pp.savefig(mpl.figure(self.label+' LID'))
        pp.close()

        mpl.figure(self.label+' ice age')
        mpl.plot(self.age_model, self.depth, color=color_mod, label='Model')
        mpl.plot(self.age, self.depth, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.age-self.sigma_age, self.age+self.sigma_age , color=color_ci)
#        mpl.plot(self.age-self.sigma_age, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_age*10, self.depth, color=color_sigma, label='$\sigma$ x10')   
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[0],self.depth[-1]))    
        mpl.legend()
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y2,y1))
        pp=PdfPages(datadir+self.label+'/ice_age.pdf')
        pp.savefig(mpl.figure(self.label+' ice age'))
        pp.close()

        mpl.figure(self.label+' gas age')
        mpl.plot(self.gage_model, self.depth, color=color_mod, label='Model')
        mpl.fill_betweenx(self.depth, self.gage-self.sigma_gage, self.gage+self.sigma_gage , color=color_ci)
        mpl.plot(self.gage, self.depth, color=color_opt, label='Corrected +/-$\sigma$')
#        mpl.plot(self.gage+self.sigma_gage, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.gage-self.sigma_gage, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_gage*10, self.depth, color=color_sigma, label='$\sigma$ x10')  
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[0],self.depth[-1]))    
        mpl.legend()
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y2,y1))
        pp=PdfPages(datadir+self.label+'/gas_age.pdf')
        pp.savefig(mpl.figure(self.label+' gas age'))
        pp.close()

        mpl.figure(self.label+' Ddepth')
        mpl.plot(self.Ddepth_model, self.depth, color=color_mod, label='Model')
        mpl.plot(self.Ddepth, self.depth, color=color_opt, label='Corrected +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.Ddepth-self.sigma_Ddepth, self.Ddepth+self.sigma_Ddepth, color=color_ci)
#        mpl.plot(self.Ddepth+self.sigma_Ddepth, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.Ddepth-self.sigma_Ddepth, self.depth, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[0],self.depth[-1]))
        mpl.legend(loc=4)
        mpl.ylim(mpl.ylim()[::-1])
        pp=PdfPages(datadir+self.label+'/Ddepth.pdf')
        pp.savefig(mpl.figure(self.label+' Ddepth'))
        pp.close()


    def save(self):
        output=np.vstack((self.depth,self.age,self.sigma_age,self.gage,self.sigma_gage,np.concatenate((self.a,np.array([self.a[-1]]))),np.concatenate((self.sigma_a,np.array([self.sigma_a[-1]]))),np.concatenate((self.tau,np.array([self.tau[-1]]))),np.concatenate((self.sigma_tau,np.array([self.sigma_tau[-1]]))),self.LID,self.sigma_LID, self.Ddepth,self.sigma_Ddepth,np.concatenate((self.a_model,np.array([self.a_model[-1]]))),np.concatenate((self.sigma_a_model,np.array([self.sigma_a_model[-1]]))),np.concatenate((self.tau_model,np.array([self.tau_model[-1]]))),np.concatenate((self.sigma_tau_model,np.array([self.sigma_tau_model[-1]]))),self.LID_model,self.sigma_LID_model))
        with open(datadir+self.label+'/output.txt','w') as f:
            f.write('#depth\tage\tsigma_age\tgas_age\tsigma_gas_age\taccu\tsigma_accu\tthinning\tsigma_thinning\tLID\tsigma_LID\tDdepth\tsigma_Ddepth\taccu_model\tsigma_accu_model\tthinning_model\tsigma_thinning_model\tLID_model\tsigma_LID_model\n')
            np.savetxt(f,np.transpose(output), delimiter='\t')
        np.savetxt(datadir+self.label+'/restart.txt',np.transpose(self.variables))
    
#    def udepth_save(self):
#        np.savetxt(datadir+self.label+'/udepth.txt',self.udepth)


class DrillingCouple:

    def __init__(self, D1, D2):
        self.D1=D1
        self.D2=D2


    def init(self):
        self.label=self.D1.label+'-'+self.D2.label
#        print 'Initialization of drilling couple ',self.label


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

        filename=datadir+self.D1.label+'-'+self.D2.label+'/gas_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.gasgasmarkers_depth1=readarray[:,0]
            self.gasgasmarkers_depth2=readarray[:,1]
            self.gasgasmarkers_sigma=readarray[:,2]
        else:
            self.gasgasmarkers_depth1=np.array([])
            self.gasgasmarkers_depth2=np.array([])
            self.gasgasmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/icegas_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.icegasmarkers_depth1=readarray[:,0]
            self.icegasmarkers_depth2=readarray[:,1]
            self.icegasmarkers_sigma=readarray[:,2]
        else:
            self.icegasmarkers_depth1=np.array([])
            self.icegasmarkers_depth2=np.array([])
            self.icegasmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/gasice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.gasicemarkers_depth1=readarray[:,0]
            self.gasicemarkers_depth2=readarray[:,1]
            self.gasicemarkers_sigma=readarray[:,2]
        else:
            self.gasicemarkers_depth1=np.array([])
            self.gasicemarkers_depth2=np.array([])
            self.gasicemarkers_sigma=np.array([])


        self.iceicemarkers_correlation=np.diag(np.ones(np.size(self.iceicemarkers_depth1)))
        self.gasgasmarkers_correlation=np.diag(np.ones(np.size(self.gasgasmarkers_depth1)))
        self.icegasmarkers_correlation=np.diag(np.ones(np.size(self.icegasmarkers_depth1)))
        self.gasicemarkers_correlation=np.diag(np.ones(np.size(self.gasicemarkers_depth1)))
        filename=datadir+'parameters-CovarianceObservations-AllDrillingCouples.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.iceicemarkers_depth1)>0:
            self.iceicemarkers_chol=np.linalg.cholesky(self.iceicemarkers_correlation)
            self.iceicemarkers_lu_piv=scipy.linalg.lu_factor(self.iceicemarkers_chol)
        if np.size(self.gasgasmarkers_depth1)>0:
            self.gasgasmarkers_chol=np.linalg.cholesky(self.gasgasmarkers_correlation)
            self.gasgasmarkers_lu_piv=scipy.linalg.lu_factor(self.gasgasmarkers_chol)
        if np.size(self.icegasmarkers_depth1)>0:
            self.icegasmarkers_chol=np.linalg.cholesky(self.icegasmarkers_correlation)
            self.icegasmarkers_lu_piv=scipy.linalg.lu_factor(self.icegasmarkers_chol)
        if np.size(self.gasicemarkers_depth1)>0:
            self.gasicemarkers_chol=np.linalg.cholesky(self.gasicemarkers_correlation)
            self.gasicemarkers_lu_piv=scipy.linalg.lu_factor(self.gasicemarkers_chol)


    def residuals(self):

        resi_iceice=(self.D1.fct_age(self.iceicemarkers_depth1)-self.D2.fct_age(self.iceicemarkers_depth2))/self.iceicemarkers_sigma
        if np.size(self.iceicemarkers_depth1)>0:
            resi_iceice=scipy.linalg.lu_solve(self.iceicemarkers_lu_piv,resi_iceice)
        resi_gasgas=(self.D1.fct_gage(self.gasgasmarkers_depth1)-self.D2.fct_gage(self.gasgasmarkers_depth2))/self.gasgasmarkers_sigma
        if np.size(self.gasgasmarkers_depth1)>0:
            resi_gasgas=scipy.linalg.lu_solve(self.gasgasmarkers_lu_piv,resi_gasgas)
        resi_icegas=(self.D1.fct_age(self.icegasmarkers_depth1)-self.D2.fct_gage(self.icegasmarkers_depth2))/self.icegasmarkers_sigma
        if np.size(self.icegasmarkers_depth1)>0:
            resi_icegas=scipy.linalg.lu_solve(self.icegasmarkers_lu_piv,resi_icegas)
        resi_gasice=(self.D1.fct_gage(self.gasicemarkers_depth1)-self.D2.fct_age(self.gasicemarkers_depth2))/self.gasicemarkers_sigma
        if np.size(self.gasicemarkers_depth1)>0:
            resi_gasice=scipy.linalg.lu_solve(self.gasicemarkers_lu_piv,resi_gasice)
        resi=np.concatenate((resi_iceice,resi_gasgas,resi_icegas,resi_gasice))
        
        return resi
    
    def display_init(self):
        
        mpl.figure(self.label+' ice-ice')
        mpl.xlabel(self.D1.label+' ice age')
        mpl.ylabel(self.D2.label+' ice age')
        if show_initial:
            mpl.errorbar(self.D1.fct_age(self.iceicemarkers_depth1),self.D2.fct_age(self.iceicemarkers_depth2), color=color_init, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")

        mpl.figure(self.label+' gas-gas')
        mpl.xlabel(self.D1.label+' gas age')
        mpl.ylabel(self.D2.label+' gas age')
        if show_initial:
            mpl.errorbar(self.D1.fct_gage(self.gasgasmarkers_depth1),self.D2.fct_gage(self.gasgasmarkers_depth2), color=color_init, xerr=self.gasgasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")

        mpl.figure(self.label+' ice-gas')
        mpl.xlabel(self.D1.label+' ice age')
        mpl.ylabel(self.D2.label+' gas age')
        if show_initial:
            mpl.errorbar(self.D1.fct_age(self.icegasmarkers_depth1),self.D2.fct_gage(self.icegasmarkers_depth2), color=color_init, xerr=self.icegasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")

        mpl.figure(self.label+' gas-ice')
        mpl.xlabel(self.D1.label+' gas age')
        mpl.ylabel(self.D2.label+' ice age')
        if show_initial:
            mpl.errorbar(self.D1.fct_gage(self.gasicemarkers_depth1),self.D2.fct_age(self.gasicemarkers_depth2), color=color_init, xerr=self.gasicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")




    def display_final(self):

        if not os.path.isdir(datadir+self.label):
            os.mkdir(datadir+self.label)


        mpl.figure(self.label+' ice-ice')
        mpl.errorbar(self.D1.fct_age_model(self.iceicemarkers_depth1),self.D2.fct_age_model(self.iceicemarkers_depth2), color=color_mod, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_age(self.iceicemarkers_depth1),self.D2.fct_age(self.iceicemarkers_depth2), color=color_opt, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc=4)
        pp=PdfPages(datadir+self.label+'/ice-ice.pdf')
        pp.savefig(mpl.figure(self.label+' ice-ice'))
        pp.close()

        mpl.figure(self.label+' gas-gas')
        mpl.errorbar(self.D1.fct_gage_model(self.gasgasmarkers_depth1),self.D2.fct_gage_model(self.gasgasmarkers_depth2), color=color_mod, xerr=self.gasgasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_gage(self.gasgasmarkers_depth1),self.D2.fct_gage(self.gasgasmarkers_depth2), color=color_opt, xerr=self.gasgasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc=4)
        pp=PdfPages(datadir+self.label+'/gas-gas.pdf')
        pp.savefig(mpl.figure(self.label+' gas-gas'))
        pp.close()

        mpl.figure(self.label+' ice-gas')
        mpl.errorbar(self.D1.fct_age_model(self.icegasmarkers_depth1),self.D2.fct_gage_model(self.icegasmarkers_depth2), color=color_mod, xerr=self.icegasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_age(self.icegasmarkers_depth1),self.D2.fct_gage(self.icegasmarkers_depth2), color=color_opt, xerr=self.icegasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc=4)
        pp=PdfPages(datadir+self.label+'/ice-gas.pdf')
        pp.savefig(mpl.figure(self.label+' ice-gas'))
        pp.close()

        mpl.figure(self.label+' gas-ice')
        mpl.errorbar(self.D1.fct_gage_model(self.gasicemarkers_depth1),self.D2.fct_age_model(self.gasicemarkers_depth2), color=color_mod, xerr=self.gasicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_gage(self.gasicemarkers_depth1),self.D2.fct_age(self.gasicemarkers_depth2), color=color_opt, xerr=self.gasicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc=4)
        pp=PdfPages(datadir+self.label+'/gas-ice.pdf')
        pp.savefig(mpl.figure(self.label+' gas-ice'))
        pp.close()

