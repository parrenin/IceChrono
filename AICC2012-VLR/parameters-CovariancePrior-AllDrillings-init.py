##In this file, the covariance matrix for the a prioris is defined.
# At first, you will most likely play with the lambda_a, lambda_LID and lambda_tau parameters which are defined parameters-AllDrillings.py file.

#Accu correlation matrix
f=interp1d(np.array([0,self.lambda_a,10000000]),np.array([1, 0, 0]))
self.correlation_corr_a=f(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age)))
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_a=gaussian(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age))/self.lambda_a)

#LID correlation matrix
f=interp1d(np.array([0,self.lambda_LID,10000000]),np.array([1, 0, 0]))
self.correlation_corr_LID=f(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age)))
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_LID=gaussian(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age))/self.lambda_LID)

#Thinning correlation matrix
g=interp1d(np.array([0,self.lambda_tau,5000]),np.array([1, 0, 0]))
self.correlation_corr_tau=g(np.abs(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth-np.transpose(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth)))
#Gaussian shape, does not work for a too high resolution
#correlation_corr_tau=gaussian(np.abs(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth-np.transpose(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth))/self.lambda_tau)


##Now we define the sigmas
# This is not used for now since the sigmas are imported from the accu-prio.txt, LID-prior.txt and thinning-prior.txt.

#def weighted_std(values, weights):
#    """
#    Return the weighted standard deviation.

#    values, weights -- Numpy ndarrays with the same shape.
#    """
#    average = np.average(values, weights=weights)
#    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
#    return (m.sqrt(variance))

#sigmaAloc=np.empty_like(self.a_model)
#mAloc=np.empty_like(self.a_model)
#for i in np.arange(np.size(self.a_model)):
#    weights=np.where((self.age_model[i]-self.age_model[:-1] <6000) & (self.age_model[i] - self.age_model[:-1]>-6000), self.age_model[1:]-self.age_model[:-1], 0)
#    sigmaAloc[i]=weighted_std(self.a_model, weights)
#    mAloc[i]=np.ma.average(self.a_model, weights=weights)

##        f=interp1d(self.depth,self.udepth_init)
#f=interp1d(self.depth,self.udepth_model, bounds_error=False, fill_value=self.udepth_model[-1]) #We should not need the bounds_error option. Check what is the problem.
#g=interp1d(self.depth[:-1], sigmaAloc, bounds_error=False, fill_value=sigmaAloc[-1])
#self.sigmap_corr_tau=self.cT1+self.cT2*f(self.corr_tau_depth) + self.cT3*g(self.corr_tau_depth)/max(sigmaAloc)

#f=interp1d(self.age_model[:-1], self.a_model, bounds_error=False, fill_value=self.a_model[-1])
#g=interp1d(self.age_model, self.depth, bounds_error=False, fill_value=self.depth[-1])
#weights=np.where(self.age_model[:-1]<12000, self.age_model[1:]-self.age_model[:-1], 0)
#A0=np.ma.average(self.a_model, weights=weights)
#self.sigmap_corr_a=self.sigmabA*np.abs(A0-f(self.corr_a_age))/max(np.abs(A0-f(self.corr_a_age)))*(1+self.cA1*g(self.corr_a_age)/self.thickness)
#self.sigmap_corr_a=np.where(self.sigmap_corr_a<self.sigmam*(1+self.cA1*g(self.corr_a_age)/self.thickness), self.sigmam*(1+self.cA1*g(self.corr_a_age)/self.thickness), self.sigmap_corr_a)

#f=interp1d(self.age_model[:-1], mAloc, bounds_error=False, fill_value=mAloc[-1])
#g=interp1d(self.corr_a_age, self.sigmap_corr_a, bounds_error=False, fill_value=self.sigmap_corr_a[-1])
#self.sigmap_corr_LID=self.sigmabL/self.sigmabA*g(self.corr_LID_age)/(1+f(self.corr_LID_age))
