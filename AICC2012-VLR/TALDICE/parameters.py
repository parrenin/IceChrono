#Parameters specific to the TALDICE ice core.
self.udepth_top=0.
self.age_top=-54.
self.depth=np.arange(0., 1486.+0.01, 1.)
self.corr_a_age=np.arange(self.age_top, 500000+self.age_top+0.01, self.age_step)
self.corr_LID_age=np.arange(self.age_top, 500000+self.age_top+0.01, self.age_step)
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))


#self.thickness=1650.
#self.cT2=0.000268
#self.sigmabA=0.6
#self.cA1=0.
#self.sigmabL=0.6

