#Parameters specific to the EDC ice core
self.udepth_top=0.                      #unthinned depth at the top of the core
self.age_top=-55.                       #age at the top of the core
self.depth=np.arange(0., 3259.3+0.01, 0.55)        #Define the grid for the age
self.thickness=3273.                    #Real thickness
self.corr_a_age=np.arange(-55, 1000000-55+0.01, self.age_step)      #Age grid for the accu correction function
self.corr_LID_age=np.arange(-55, 1000000-55+0.01, self.age_step)    #Age grid for the LID correction function
self.cT2=0.000030
self.sigmabA=0.7
self.cA1=0.
self.sigmabL=0.7
self.restart=False
