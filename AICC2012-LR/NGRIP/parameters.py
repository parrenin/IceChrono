#Parameters specific to the NGRIP ice core
self.udepth_top=3.313966267738
self.age_top=-30.
self.depth=np.arange(8., 3084.+0.01, 1.)
self.thickness=3085. #From NGRIP community members (2004)
self.corr_a_age=np.arange(-30, 150000-30+0.01, self.age_step)
self.corr_LID_age=np.arange(-30, 150000-30+0.01, self.age_step)
self.cT2=0.000016
self.sigmabA=0.6
self.cA1=0.
self.sigmabL=0.6
self.restart=False
