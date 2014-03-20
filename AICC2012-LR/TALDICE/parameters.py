#Parameters specific to the TALDICE ice core, from the Gnumeric TALDICE file
#Changelog: modified the udepth.txt file to prevent an error during the leastsq function
#TODO: use the real EDC-TALDICE ice-ice strati links with small errors by using a restart file
#TODO: TALDICE does not run alone. Check what is the problem.

self.dim=1
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.calc_udepth=True
self.udepth_min=0.
self.age_min=-54.
self.depth_min=0.
self.depth_max=1486.
self.step=1.
self.gamma_source=3.4 #Same as DC, need to check with B. Stenni if there is a better value
self.beta_source=1.5  #Same as DC, need to check with B. Stenni if there is a better value
self.A0=8.26e-02
self.beta=6.64e-03
self.thickness=1620. #1650 is the documented value but we need 1620 for the code to not crash
self.pprime=9.06e-01
self.s=-1.81e-01
self.mu=1.57e-02
self.age_max=500000.+self.age_min
self.corr_a=np.zeros(51)
self.corr_LID=np.zeros(51)
self.LID_value=76.
self.Dfirn=0.698    #This is the EDC value
self.cT2=0.000268
self.sigmabA=0.6
self.cA1=0.
self.restart=False
