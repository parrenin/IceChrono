##In this file, the covariance matrix for the a prioris is defined.
# At first, you will most likely play with the lambda_a, lambda_LID and lambda_tau parameters which are defined parameters-AllDrillings.py file.

#Accu correlation matrix
self.correlation_corr_a=np.interp(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age)), np.array([0,self.lambda_a]),np.array([1, 0]))
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_a=gaussian(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age))/self.lambda_a)

#LID correlation matrix
self.correlation_corr_LID=np.interp(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age)), np.array([0,self.lambda_LID]),np.array([1, 0]) )
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_LID=gaussian(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age))/self.lambda_LID)

#Thinning correlation matrix
self.correlation_corr_tau=np.interp(np.abs(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth-np.transpose(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth)), np.array([0,self.lambda_tau]),np.array([1, 0]) )
#Gaussian shape, does not work for a too high resolution
#correlation_corr_tau=gaussian(np.abs(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth-np.transpose(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth))/self.lambda_tau)

