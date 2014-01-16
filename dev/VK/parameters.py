#Parameters specific to the Vostok ice core
self.dim=2
self.calc_a=False
self.calc_tau=False
self.calc_LID=True
self.age_min=-47
self.min_depth=0.
self.max_depth=3501.
self.thickness=3767.
self.step=1
#self.gamma_source=3.4
#self.beta_source=1.5
self.lambda_a=4000
self.sigmap_corr_a=0.2
self.k=0.1
self.lambda_tau=70
self.lambda_LIDIE=4000
self.sigmap_corr_LIDIE=0.4
self.age_max=800000.+self.age_min
self.corr_a=np.zeros(81)
self.corr_LIDIE=np.zeros(81)
#self.LID_depth=np.array([0., (self.max_depth-0.1)])
#self.LID_LID=np.array([98., 98.])  #TODO: change this EDC value
self.LID_value=98.   #TODO: change this EDC value
self.corr_tau=np.zeros(101)
self.restart=False


#self.pprime=1.19
#self.s=0.
#self.mu=0.00066/self.A0
#self.beta=0.0157
#self.A0=0.02841
