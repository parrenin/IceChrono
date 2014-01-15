#Parameters specific to the EDC ice core
#TODO: Check these parameters with Parrenin et al. (CP, 2007), especially ice thickness
self.dim=1
self.calc_a=True
self.calc_tau=True
self.age_surf=-50.
self.max_depth=3190.
self.step=0.55
self.gamma_source=3.4
self.beta_source=1.5
self.A0=3.30e-02
self.lambda_a=4000
self.beta=1.65e-02
self.sigmap_corr_a=0.2
self.thickness=3273.
self.pprime=1.59
self.s=-2.82e-01
self.mu=5.34e-02
self.k=0.1
self.lambda_tau=70
self.lambda_LIDIE=4000
self.sigmap_corr_LIDIE=0.4
self.age_max=1000000.
self.corr_a=np.zeros(101)
self.corr_LIDIE=np.zeros(101)
self.LID_depth=np.array([0., self.max_depth])
self.LID_LID=np.array([98., 98.])
self.corr_tau=np.zeros(101)
self.restart=False


