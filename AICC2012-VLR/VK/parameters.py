#Parameters specific to the Vostok ice core
self.udepth_top=0.
self.age_top=-50.
self.depth=np.arange(0., 3501+0.01, 1.)
self.corr_a_age=np.arange(self.age_top, 800000+self.age_top+0.01, self.age_step)
self.corr_LID_age=np.arange(self.age_top, 800000+self.age_top+0.01, self.age_step)
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))


#self.thickness=3767.
#self.cT2=0.000084
#self.sigmabA=0.6
#self.cA1=1.
#self.sigmabL=0.7

