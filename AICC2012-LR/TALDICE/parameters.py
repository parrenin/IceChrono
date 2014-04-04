#Parameters specific to the TALDICE ice core.
self.udepth_min=0.
self.age_min=-54.
self.depth_min=0.
self.depth_max=1486.
self.step=1.
self.thickness=1650.
self.age_max=500000.+self.age_min
self.corr_a=np.zeros(51)
self.corr_LID=np.zeros(51)
self.Dfirn=0.698    #This is the EDC value
self.cT2=0.000268
self.sigmabA=0.6
self.cA1=0.
self.sigmabL=0.6
self.restart=False
