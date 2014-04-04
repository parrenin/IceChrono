#Parameters specific to the EDML ice core
self.udepth_top=8.54025
self.age_top=70.
self.depth_top=18.
self.depth_bot=2564.
self.thickness=3000
self.step=1
self.depth=np.arange(self.depth_top, self.depth_bot+0.01, self.step)
self.age_bot=300000.+self.age_top
self.corr_a=np.zeros(31)
self.corr_LID=np.zeros(31)
self.cT2=0.000078
self.sigmabA=0.5
self.cA1=1.
self.sigmabL=0.6
self.restart=False

