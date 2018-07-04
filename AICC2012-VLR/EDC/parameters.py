#Parameters specific to the EDC ice core
self.udepth_top=0.                      #unthinned depth at the top of the core
self.age_top=-55.                       #age at the top of the core
self.depth=np.arange(0., 3259.3+0.01, 0.55)        #Define the depth grid for the age calculation
self.corr_a_age=np.arange(self.age_top, 1000000+self.age_top+0.01, self.age_step)      #Age grid for the accu correction function
self.corr_LID_age=np.arange(self.age_top, 1000000+self.age_top+0.01, self.age_step)    #Age grid for the LID correction function
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))   #Depth grid for the thinning correction function
self.accu_prior_rep='staircase'    #linear or staircase. Define whether the prior accu representation is linear or staircase in-between the data points.


#The following parameters defines the covariance matrices as in AICC2012 (Bazin et al., 2013 and Veres et al., 2013).
#self.thickness=3273.                    #Real thickness
#self.cT2=0.000030/0.55
#self.sigmabA=0.7
#self.cA1=0.
#self.sigmabL=0.7
