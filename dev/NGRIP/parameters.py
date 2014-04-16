#Parameters specific to the NGRIP ice core
self.dim=1
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.udepth_top=3.313966267738
self.age_top=-30.
self.depth=np.arange(8., 3084.+0.01, 1.)
self.thickness=3085. #From NGRIP community members (2004)
self.corr_a_age=np.arange(-30, 150000-30+0.01, self.age_step)
self.corr_LID_age=np.arange(-30, 150000-30+0.01, self.age_step)
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))
self.restart=False
