#Parameters specific to the EDML ice core
#TODO: udepth.txt is incorrect
self.dim=2
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.calc_udepth_init=False
self.age_min=-47
self.depth_min=0.
self.depth_max=2416.
self.thickness=3000
self.step=1
#self.gamma_source=3.4
#self.beta_source=1.5
self.lambda_a=4000
self.sigmap_corr_a=0.2
self.k=0.1
self.lambda_tau=70
self.lambda_LID=4000
self.sigmap_corr_LID=0.4
self.age_max=300000.+self.age_min
self.corr_a=np.zeros(31)
self.corr_LID=np.zeros(31)
self.LID_value=98.
self.Dfirn=0.698    #This is the EDC value
self.corr_tau=np.zeros(101)
self.restart=False


#self.thickness=3237
#self.pprime=1.19
#self.s=0.
#self.mu=0.00066/self.A0
#self.beta=0.0157
#self.A0=0.02841
