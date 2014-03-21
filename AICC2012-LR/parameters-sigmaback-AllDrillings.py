def weighted_std(values, weights):
    """
    Return the weighted standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (m.sqrt(variance))

sigmaAloc=np.empty_like(self.a_model)
mAloc=np.empty_like(self.a_model)
for i in np.arange(np.size(self.a_model)):
    weights=np.where((self.age_model[i]-self.age_model[:-1] <6000) & (self.age_model[i] - self.age_model[:-1]>-6000), self.age_model[1:]-self.age_model[:-1], 0)
    sigmaAloc[i]=weighted_std(self.a_model, weights)
    mAloc[i]=np.ma.average(self.a_model, weights=weights)

#        f=interpolate.interp1d(self.depth,self.udepth_init)
f=interpolate.interp1d(self.depth,self.udepth_init, bounds_error=False, fill_value=self.udepth_init[-1]) #We should not need the bounds_error option. Check what is the problem.
g=interpolate.interp1d(self.depth[:-1], sigmaAloc, bounds_error=False, fill_value=sigmaAloc[-1])
self.sigmap_corr_tau=self.cT1+self.cT2*f(self.corr_tau_depth) + self.cT3*g(self.corr_tau_depth)/max(sigmaAloc)

f=interpolate.interp1d(self.age_model[:-1], self.a_model, bounds_error=False, fill_value=self.a_model[-1])
g=interpolate.interp1d(self.age_model, self.depth, bounds_error=False, fill_value=self.depth[-1])
weights=np.where(self.age_model[:-1]<12000, self.age_model[1:]-self.age_model[:-1], 0)
A0=np.ma.average(self.a_model, weights=weights)
self.sigmap_corr_a=self.sigmabA*np.abs(A0-f(self.corr_a_age))/max(np.abs(A0-f(self.corr_a_age)))*(1+self.cA1*g(self.corr_a_age)/self.thickness)
self.sigmap_corr_a=np.where(self.sigmap_corr_a<self.sigmam, self.sigmam*(1+self.cA1*g(self.corr_a_age)/self.thickness), self.sigmap_corr_a)

f=interpolate.interp1d(self.age_model[:-1], mAloc, bounds_error=False, fill_value=mAloc[-1])
g=interpolate.interp1d(self.corr_a_age, self.sigmap_corr_a, bounds_error=False, fill_value=self.sigmap_corr_a[-1])
self.sigmap_corr_LID=self.sigmabL/self.sigmabA*g(self.corr_LID_age)/(1+f(self.corr_LID_age)/max(mAloc))