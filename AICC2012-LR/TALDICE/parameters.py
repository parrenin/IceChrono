#Parameters specific to the TALDICE ice core.
self.udepth_top=0.
self.age_top=-54.
self.depth_top=0.
self.depth_bot=1486.
self.step=1.
self.depth=np.arange(self.depth_top, self.depth_bot+0.01, self.step)
self.thickness=1650.
self.age_bot=500000.+self.age_top
self.corr_a=np.zeros(51)
self.corr_LID=np.zeros(51)
self.cT2=0.000268
self.sigmabA=0.6
self.cA1=0.
self.sigmabL=0.6
self.restart=False
