#Parameters specific to the Vostok ice core
self.dim=2
self.udepth_min=0.
self.age_min=-50.
self.depth_min=0.
self.depth_max=3501.
self.thickness=3767.
self.step=1
#self.gamma_source=3.4
#self.beta_source=1.5
self.age_max=800000.+self.age_min
self.corr_a=np.zeros(81)
self.corr_LID=np.zeros(81)
self.LID_value=98.   #TODO: change this EDC value
self.Dfirn=0.698    #This is the EDC value
self.cT2=0.000084
self.sigmabA=0.6
self.cA1=1.
self.sigmabL=0.7
self.restart=False


#self.pprime=1.19
#self.s=0.
#self.mu=0.00066/self.A0
#self.beta=0.0157
#self.A0=0.02841
