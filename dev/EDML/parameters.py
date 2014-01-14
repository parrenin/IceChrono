#Parameters specific to the EDML ice core
#TODO: tau_model.txt is incorrect, takes into account firn compaction
#TODO: udepth.txt is incorrect
self.dim=2
self.calc_a=False
self.calc_tau=False
self.age_surf=-47
self.max_depth=2416.
self.thickness=3000
self.step=1
#self.gamma_source=3.4
#self.beta_source=1.5
self.lambda_a=4000
self.sigmap_corr_a=0.2
self.k=0.1
self.lambda_tau=70
self.lambda_LIDIE=4000
self.sigmap_corr_LIDIE=0.4
self.age_max=800000.
self.corr_a=np.zeros(31)
self.corr_LIDIE=np.zeros(31)
self.LID_depth=np.array([0., (self.max_depth-0.1)])
self.LID_LID=np.array([98., 98.])
self.corr_tau=np.zeros(101)
self.restart=True


#self.thickness=3237
#self.pprime=1.19
#self.s=0.
#self.mu=0.00066/self.A0
#self.beta=0.0157
#self.A0=0.02841
