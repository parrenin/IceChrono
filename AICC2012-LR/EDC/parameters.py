#Parameters specific to the EDC ice core
self.udepth_min=0.
self.age_min=-55.
self.depth_min=0.
self.depth_max=3259.3
self.step=0.55
self.thickness=3273.
self.age_max=1000000.+self.age_min
self.corr_a=np.zeros(101)
self.corr_LID=np.zeros(101)
self.Dfirn=0.698    #From Parrenin et al. (CP, 2012b)
self.cT2=0.000030
self.sigmabA=0.7
self.cA1=0.
self.sigmabL=0.7
self.restart=False
