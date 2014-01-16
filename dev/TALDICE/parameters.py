#Parameters specific to the TALDICE ice core, from the Gnumeric TALDICE file
#Changelog: modified the udepth.txt file to prevent an error during the leastsq function
#TODO: use the real EDC-TALDICE ice-ice strati links with small errors by using a restart file
#TODO: TALDICE does not run alone. Check what is the problem.

self.dim=1
self.calc_a=True
self.calc_tau=True
self.age_min=-50
self.min_depth=0.
self.max_depth=1620.
self.step=1.
self.gamma_source=3.4 #Same as DC, need to check with B. Stenni if there is a better value
self.beta_source=1.5  #Same as DC, need to check with B. Stenni if there is a better value
self.A0=8.26e-02
self.lambda_a=4000
self.beta=6.64e-03
self.sigmap_corr_a=0.2
self.thickness=1620. #1650 is the documented value but we need 1620 for the code to not crash
self.pprime=9.06e-01
self.s=-1.81e-01
self.mu=1.57e-02
self.k=0.1
self.lambda_tau=70
self.lambda_LIDIE=4000
self.sigmap_corr_LIDIE=0.4
self.age_max=300000.
self.corr_a=np.zeros(51)
self.corr_LIDIE=np.zeros(51)
self.LID_depth=np.array([0., self.max_depth])
self.LID_LID=np.array([76., 76.])
self.corr_tau=np.zeros(101)
self.restart=False
