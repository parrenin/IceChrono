#Parameters specific to the EDML ice core
self.udepth_top=8.54025
self.age_top=70.
self.depth=np.arange(18., 2564.+0.01, 1.)
self.corr_a_age=np.arange(70., 300000+70.+0.01, self.age_step)
self.corr_LID_age=np.arange(70., 300000+70.+0.01, self.age_step)
self.thickness=3000
self.cT2=0.000078
self.sigmabA=0.5
self.cA1=1.
self.sigmabL=0.6
self.restart=False

