#        f=interpolate.interp1d(self.depth,self.udepth_init)
f=interpolate.interp1d(self.depth,self.udepth_init, bounds_error=False, fill_value=self.udepth_init[-1]) #We should not need the bounds_error option. Check what is the problem.
self.sigmap_corr_tau=self.cT1+self.cT2*f(self.corr_tau_depth)

#f=interpolate.interp1d(self.age_model[:-1], self.a_model, bounds_error=False, fill_value=self.a_model[-1])
#g=interpolate.interp1d(self.age_model, self.depth, bounds_error=False, fill_value=self.depth[-1])
#self.sigmap_corr_a=self.sigmabA*np.abs(self.A0-f(self.corr_a_age))/max(np.abs(self.A0-f(self.corr_a_age)))*(1+self.cA1*g(self.corr_a_age)/self.thickness)
#self.sigmap_corr_a=np.where(self.sigmap_corr_a<self.sigmam, self.sigmam*(1+self.cA1*g(self.corr_a_age)/self.thickness), self.sigmap_corr_a)
