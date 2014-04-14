#Parameters specific to the TALDICE ice core.
self.udepth_top=0.
self.age_top=-54.
self.depth=np.arange(0., 1486.+0.01, 1.)
self.corr_a_age=np.arange(-50, 500000-50+0.01, self.age_step)
self.corr_LID_age=np.arange(-50, 500000-50+0.01, self.age_step)
self.thickness=1650.
self.cT2=0.000268
self.sigmabA=0.6
self.cA1=0.
self.sigmabL=0.6
self.restart=False
