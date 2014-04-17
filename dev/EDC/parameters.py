#Parameters specific to the EDC ice core
#TODO: Check these parameters with Parrenin et al. (CP, 2007), especially ice thickness
self.dim=1
self.udepth_top=0.
self.calc_a=True
self.calc_a_method='fullcorr'
self.calc_tau=True
self.calc_LID=True
self.age_top=-50.
self.depth=np.arange(0., 3190.+0.01, 0.55)
self.gamma_source=3.4
self.beta_source=1.5
self.A0=3.30e-02
self.beta=1.65e-02
self.thickness=3273.
self.pprime=m.log(2.3+1)
self.s=0.
self.muprime=m.log(5.34e-02)
self.corr_a_age=np.arange(-55, 1000000-55+0.01, self.age_step)      #Age grid for the accu correction function
self.corr_LID_age=np.arange(-55, 1000000-55+0.01, self.age_step)    #Age grid for the LID correction function
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))
self.LID_value=98.
self.restart=False
