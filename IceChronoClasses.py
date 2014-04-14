#TODO: extend the chronology down to the bedrock by extrapolating the accumulation
#TODO: optinally use a restart file to have a bootstrap method
#TODO: is there an elegant way to unpack the variables vector in the model function?
#TODO: allow to save the correction vector to be able to restart while changing the resolution
#TODO: include some checks for when dDdepth/dz>1


def gaussian(x):
    return np.exp(-x**2/2)


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
            self.depthtop=readarray[:,0]
            self.d18Oice=readarray[:,1]
            self.deutice=readarray[:,2]
            self.d18Osw=readarray[:,3]
            self.excess=self.deutice-8*self.d18Oice   # dans Uemura : d=excess
            self.a=np.empty_like(self.deutice)
            self.d18Oice_corr=self.d18Oice-self.d18Osw*(1+self.d18Oice/1000)/(1+self.d18Osw/1000)	#Uemura (1)
            self.deutice_corr=self.deutice-8*self.d18Osw*(1+self.deutice/1000)/(1+8*self.d18Osw/1000)	#Uemura et al. (CP, 2012) (2) 
            self.excess_corr=self.deutice_corr-8*self.d18Oice_corr
            self.deutice_fullcorr=self.deutice_corr+self.gamma_source/self.beta_source*self.excess_corr
        else:
            self.a_model=np.loadtxt(datadir+self.label+'/accu-prior.txt')
            self.a=self.a_model

        
        self.age=np.empty_like(self.depth)
        self.gage=np.empty_like(self.depth)
        

        readarray=np.loadtxt(datadir+self.label+'/density-prior.txt')
#        self.density_depth=readarray[:,0]
        self.D=readarray
        self.iedepth=np.cumsum(np.concatenate((np.array([0]), self.D*self.depth_inter)))
        self.thickness_ie=self.thickness-self.depth[-1]+self.iedepth[-1]
        
        if self.calc_LID:
            if self.depth[0]<self.LID_value:
                self.LID_depth=np.array([self.depth[0], self.LID_value, self.depth[-1]])
                self.LID_LID=np.array([self.depth[0], self.LID_value, self.LID_value])
            else:
                self.LID_depth=np.array([self.depth[0], self.depth[-1]])
                self.LID_LID=np.array([self.LID_value, self.LID_value])
            f=interpolate.interp1d(self.LID_depth, self.LID_LID, bounds_error=False, fill_value=self.LID_LID[-1])
            self.LID_model=f(self.depth)
        else:
            self.LID_model=np.loadtxt(datadir+self.label+'/LID-prior.txt')



        self.Ddepth=np.empty_like(self.depth)
        self.udepth=np.empty_like(self.depth)

#        print 'depth_mid ', np.size(self.depth_mid)
        self.zeta=(self.thickness-self.depth_mid)/self.thickness
#        print 'zeta ', np.size(self.zeta)
        if self.calc_tau:
            self.tau=np.empty_like(self.depth_mid)
        else:
            self.tau_model=np.loadtxt(datadir+self.label+'/thinning-prior.txt')
            self.tau=self.tau_model

        self.raw_model()

## Now we set up the correction functions


        self.corr_a_age=np.arange(self.age_top, self.age_bot+0.1, (self.age_bot-self.age_top)/(np.size(self.corr_a)-1))
        self.correlation_corr_a=np.empty((np.size(self.corr_a),np.size(self.corr_a)))



        self.corr_LID=np.zeros(np.size(self.corr_LID))
        self.corr_LID_age=np.arange(self.age_top,self.age_bot+0.1, (self.age_bot-self.age_top)/(np.size(self.corr_LID)-1))
        self.correlation_corr_LID=np.empty((np.size(self.corr_LID),np.size(self.corr_LID)))


        self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(np.size(self.corr_tau)-1))
        self.corr_tau=np.zeros(np.size(self.corr_tau))
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

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=np.linalg.cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=np.linalg.cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=np.linalg.cholesky(self.correlation_corr_tau)


        if not self.calc_a_sigma:
            readarray=np.loadtxt(datadir+self.label+'/accu-sigma-prior.txt')
            f=interpolate.interp1d(self.age_model[:-1],readarray, bounds_error=False, fill_value=readarray[-1])
            self.sigmap_corr_a=f(self.corr_a_age)
        if not self.calc_tau_sigma:
            readarray=np.loadtxt(datadir+self.label+'/thinning-sigma-prior.txt')
            f=interpolate.interp1d(self.depth[:-1],readarray, bounds_error=False, fill_value=readarray[-1])
            self.sigmap_corr_tau=f(self.corr_tau_depth)
        if not self.calc_LID_sigma:
            readarray=np.loadtxt(datadir+self.label+'/LID-sigma-prior.txt')
            f=interpolate.interp1d(self.age_model,readarray, bounds_error=False, fill_value=readarray[-1])
            self.sigmap_corr_LID=f(self.corr_LID_age)


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
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.icemarkers_depth=readarray[:,0]
            self.icemarkers_age=readarray[:,1]
            self.icemarkers_sigma=readarray[:,2]
        else:
            self.icemarkers_depth=np.array([])
            self.icemarkers_age=np.array([])
            self.icemarkers_sigma=np.array([])

        filename=datadir+self.label+'/gas_age.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.gasmarkers_depth=readarray[:,0]
            self.gasmarkers_age=readarray[:,1]
            self.gasmarkers_sigma=readarray[:,2]
        else:
            self.gasmarkers_depth=np.array([])
            self.gasmarkers_age=np.array([])
            self.gasmarkers_sigma=np.array([])

        filename=datadir+self.label+'/ice_age_intervals.txt'
        if os.path.isfile(filename) and open(filename).read():
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
        if os.path.isfile(filename) and open(filename).read():
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
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.Ddepth_depth=readarray[:,0]
            self.Ddepth_Ddepth=readarray[:,1]
            self.Ddepth_sigma=readarray[:,2]
        else:
            self.Ddepth_depth=np.array([])
            self.Ddepth_Ddepth=np.array([])
            self.Ddepth_sigma=np.array([])


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
        
        g_model=interpolate.interp1d(self.iedepth, self.udepth_model)
        self.LIDIE_model=self.LID_model*self.Dfirn
        self.ULIDIE_model=g_model(self.LIDIE_model)
        i_model=interpolate.interp1d(self.udepth_model, self.depth)

        #Ice age
        self.icelayerthick_model=self.tau_model*self.a_model/self.D
        self.age_model=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model/self.a_model*self.depth_inter)))
            
        f_model=interpolate.interp1d(self.depth, self.age_model, bounds_error=False, fill_value=np.nan)

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

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=np.linalg.cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=np.linalg.cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=np.linalg.cholesky(self.correlation_corr_tau)


        #Accu
        corr=np.dot(self.chol_a,self.corr_a)*self.sigmap_corr_a
        j=interpolate.interp1d(self.corr_a_age, corr, bounds_error=False, fill_value=corr[-1])
        self.a=self.a_model*np.exp(j(self.age_model[:-1])) #FIXME: we should use mid-age and not age

        #Thinning
        h=interpolate.interp1d(self.corr_tau_depth, np.dot(self.chol_tau,self.corr_tau)*self.sigmap_corr_tau)
        self.tau=self.tau_model*np.exp(h(self.depth_mid))
        self.udepth=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau*self.depth_inter)))
        g=interpolate.interp1d(self.iedepth, self.udepth)
        corr=np.dot(self.chol_LID,self.corr_LID)*self.sigmap_corr_LID
        j=interpolate.interp1d(self.corr_LID_age, corr, bounds_error=False, fill_value=corr[-1])
        self.LID=self.LID_model*np.exp(j(self.age_model))
        self.LIDIE=self.LID*self.Dfirn
        self.ULIDIE=g(self.LIDIE)
        i=interpolate.interp1d(self.udepth, self.depth)

        #Ice age
        self.icelayerthick=self.tau*self.a/self.D
        self.age=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau/self.a*self.depth_inter)))
        f=interpolate.interp1d(self.depth,self.age, bounds_error=False, fill_value=np.nan)

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
        f=interpolate.interp1d(self.depth,self.age)
        return f(depth)
   
    def fct_age_model(self, depth):
        f=interpolate.interp1d(self.depth,self.age_model)
        return f(depth)
   
    def fct_gage(self, depth):
        f=interpolate.interp1d(self.depth,self.gage)
        return f(depth)

    def fct_gage_model(self, depth):
        f=interpolate.interp1d(self.depth,self.gage_model)
        return f(depth)

    def fct_Ddepth(self, depth):
        f=interpolate.interp1d(self.depth,self.Ddepth)
        return f(depth)

    def residuals(self, variables):
        self.model(variables)
        resi_age=(self.fct_age(self.icemarkers_depth)-self.icemarkers_age)/self.icemarkers_sigma
        resi_gage=(self.fct_gage(self.gasmarkers_depth)-self.gasmarkers_age)/self.gasmarkers_sigma
        resi_Ddepth=(self.fct_Ddepth(self.Ddepth_depth)-self.Ddepth_Ddepth)/self.Ddepth_sigma
        resi_iceint=(self.fct_age(self.iceintervals_depthbot)-self.fct_age(self.iceintervals_depthtop)-self.iceintervals_duration)/self.iceintervals_sigma
        resi_gasint=(self.fct_gage(self.gasintervals_depthbot)-self.fct_gage(self.gasintervals_depthtop)-self.gasintervals_duration)/self.gasintervals_sigma
        resi_corr_tau=self.corr_tau
        resi_corr_a=self.corr_a
        resi_corr_LID=self.corr_LID
        return np.concatenate((resi_age,resi_gage, resi_Ddepth, resi_iceint, resi_gasint, resi_corr_tau, resi_corr_a, resi_corr_LID))


    def jacobian(self):
        epsilon=np.sqrt(np.diag(self.hess))/100000000.
        model0=self.model(self.variables)
        jacob=np.empty((np.size(self.variables), np.size(model0)))
        for i in np.arange(np.size(self.variables)):
            var=self.variables+0
            var[i]=var[i]+epsilon[i]
            model1=self.model(var)
            jacob[i]=(model1-model0)/epsilon[i]
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
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.age)]),np.dot(self.hess,jacob[:,index:index+np.size(self.age)]))
        self.sigma_age=np.sqrt(np.diag(c_model))
        index=index+np.size(self.age)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.gage)]),np.dot(self.hess,jacob[:,index:index+np.size(self.gage)]))
        self.sigma_gage=np.sqrt(np.diag(c_model))
        index=index+np.size(self.gage)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.Ddepth)]),np.dot(self.hess,jacob[:,index:index+np.size(self.Ddepth)]))
        self.sigma_Ddepth=np.sqrt(np.diag(c_model))
        index=index+np.size(self.Ddepth)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.a)]),np.dot(self.hess,jacob[:,index:index+np.size(self.a)]))
        self.sigma_a=np.sqrt(np.diag(c_model))
        index=index+np.size(self.a)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.tau)]),np.dot(self.hess,jacob[:,index:index+np.size(self.tau)]))
        self.sigma_tau=np.sqrt(np.diag(c_model))
        index=index+np.size(self.tau)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.LID)]),np.dot(self.hess,jacob[:,index:index+np.size(self.LID)]))
        self.sigma_LID=np.sqrt(np.diag(c_model))
        index=index+np.size(self.LID)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.icelayerthick)]),np.dot(self.hess,jacob[:,index:index+np.size(self.icelayerthick)]))
        self.sigma_icelayerthick=np.sqrt(np.diag(c_model))
        index=index+np.size(self.icelayerthick)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.gaslayerthick)]),np.dot(self.hess,jacob[:,index:index+np.size(self.gaslayerthick)]))
        self.sigma_gaslayerthick=np.sqrt(np.diag(c_model))
        index=index+np.size(self.gaslayerthick)
        
        return self.sigma_age, self.sigma_gage, self.sigma_Ddepth, self.sigma_a, self.sigma_tau, self.sigma_LID, self.sigma_icelayerthick

    

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
        output=np.vstack((self.depth,self.age,self.sigma_age,self.gage,self.sigma_gage,np.concatenate((self.a,np.array([self.a[-1]]))),np.concatenate((self.sigma_a,np.array([self.sigma_a[-1]]))),np.concatenate((self.tau,np.array([self.tau[-1]]))),np.concatenate((self.sigma_tau,np.array([self.sigma_tau[-1]]))),self.LID,self.sigma_LID, self.Ddepth,self.sigma_Ddepth))
        np.savetxt(datadir+self.label+'/output.txt',np.transpose(output), header='depth age sigma_age gas_age sigma_gas_age accu sigma_accu thinning sigma_thinning LID sigma_LID Ddepth sigma_Ddepth')
        np.savetxt(datadir+self.label+'/restart.txt',np.transpose(self.variables))
    
    def udepth_save(self):
        np.savetxt(datadir+self.label+'/udepth.txt',self.udepth)
        return

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
            self.icemarkers_depth1=readarray[:,0]
            self.icemarkers_depth2=readarray[:,1]
            self.icemarkers_sigma=readarray[:,2]
        else:
            self.icemarkers_depth1=np.array([])
            self.icemarkers_depth2=np.array([])
            self.icemarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/gas_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.gasmarkers_depth1=readarray[:,0]
            self.gasmarkers_depth2=readarray[:,1]
            self.gasmarkers_sigma=readarray[:,2]
        else:
            self.gasmarkers_depth1=np.array([])
            self.gasmarkers_depth2=np.array([])
            self.gasmarkers_sigma=np.array([])

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

        return

    def residuals(self):

        resi_ice=(self.D1.fct_age(self.icemarkers_depth1)-self.D2.fct_age(self.icemarkers_depth2))/self.icemarkers_sigma
        resi_gas=(self.D1.fct_gage(self.gasmarkers_depth1)-self.D2.fct_gage(self.gasmarkers_depth2))/self.gasmarkers_sigma
        resi_icegas=(self.D1.fct_age(self.icegasmarkers_depth1)-self.D2.fct_gage(self.icegasmarkers_depth2))/self.icegasmarkers_sigma
        resi_gasice=(self.D1.fct_gage(self.gasicemarkers_depth1)-self.D2.fct_age(self.gasicemarkers_depth2))/self.gasicemarkers_sigma
        resi=np.concatenate((resi_ice,resi_gas,resi_icegas,resi_gasice))
        
        return resi
    
    def display_init(self):
        
        mpl.figure(self.label+' ice-ice')
        mpl.xlabel(self.D1.label+' ice age')
        mpl.ylabel(self.D2.label+' ice age')
        if show_initial:
            mpl.errorbar(self.D1.fct_age(self.icemarkers_depth1),self.D2.fct_age(self.icemarkers_depth2), color=color_init, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")

        mpl.figure(self.label+' gas-gas')
        mpl.xlabel(self.D1.label+' gas age')
        mpl.ylabel(self.D2.label+' gas age')
        if show_initial:
            mpl.errorbar(self.D1.fct_gage(self.gasmarkers_depth1),self.D2.fct_gage(self.gasmarkers_depth2), color=color_init, xerr=self.gasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")

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

        return



    def display_final(self):

        if not os.path.isdir(datadir+self.label):
            os.mkdir(datadir+self.label)


        mpl.figure(self.label+' ice-ice')
        mpl.errorbar(self.D1.fct_age_model(self.icemarkers_depth1),self.D2.fct_age_model(self.icemarkers_depth2), color=color_mod, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_age(self.icemarkers_depth1),self.D2.fct_age(self.icemarkers_depth2), color=color_opt, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
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
        mpl.errorbar(self.D1.fct_gage_model(self.gasmarkers_depth1),self.D2.fct_gage_model(self.gasmarkers_depth2), color=color_mod, xerr=self.gasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Model")
        mpl.errorbar(self.D1.fct_gage(self.gasmarkers_depth1),self.D2.fct_gage(self.gasmarkers_depth2), color=color_opt, xerr=self.gasmarkers_sigma, linestyle='', marker='o', markersize=2, label="Optimized")
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
        
        return
