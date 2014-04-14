#Parameters specific to the EDC ice core
self.udepth_top=0.                      #unthinned depth at the top of the core
self.age_top=-55.                       #age at the top of the core
self.depth=np.arange(0., 3259.3+0.01, 0.55)        #Define the grid for the age
self.thickness=3273.                    #Real thickness
self.age_bot=1000000.+self.age_top      #Upper bound of the age at bottom (for the definition of the accu and thinning correction functions)
self.corr_a=np.zeros(101)               #Number of nodes for accu correction
self.corr_LID=np.zeros(101)             #Number of nodes for LID correction
self.cT2=0.000030
self.sigmabA=0.7
self.cA1=0.
self.sigmabL=0.7
self.restart=False
