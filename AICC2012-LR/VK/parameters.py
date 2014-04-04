#Parameters specific to the Vostok ice core
self.udepth_top=0.
self.age_top=-50.
self.depth_top=0.
self.depth_bot=3501.
self.thickness=3767.
self.step=1
self.depth=np.arange(self.depth_top, self.depth_bot+0.01, self.step)
self.age_bot=800000.+self.age_top
self.corr_a=np.zeros(81)
self.corr_LID=np.zeros(81)
self.cT2=0.000084
self.sigmabA=0.6
self.cA1=1.
self.sigmabL=0.7
self.restart=False

