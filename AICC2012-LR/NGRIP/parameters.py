#Parameters specific to the NGRIP ice core
self.udepth_top=3.313966267738
self.age_top=-30.
self.depth_top=8.
self.depth_bot=3084.
self.step=1
self.thickness=3085. #From NGRIP community members (2004)
self.age_bot=150000.+self.age_top
self.corr_a=np.zeros(16)
self.corr_LID=np.zeros(16)
self.cT2=0.000016
self.sigmabA=0.6
self.cA1=0.
self.sigmabL=0.6
self.restart=False
