#Parameters specific to the EDML ice core
self.udepth_top=8.54025
self.age_top=70.
self.depth=np.arange(18., 2564.+0.01, 1.)
self.age_bot=300000.+self.age_top
self.corr_a=np.zeros(31)
self.corr_LID=np.zeros(31)
self.thickness=3000
self.cT2=0.000078
self.sigmabA=0.5
self.cA1=1.
self.sigmabL=0.6
self.restart=False

