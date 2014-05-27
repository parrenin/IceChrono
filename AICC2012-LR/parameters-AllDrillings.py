self.calc_a=False           #Use False for now.
self.calc_tau=False         #Use False for now.
self.calc_LID=False         #Use False for now.
self.corr_tau_nodes=51  #Define the number of nodes of the thinning function.
self.age_step=10000.	#Define the age step for the LID and accu correction functions.
self.Dfirn=0.7            #Average density of the firn

#Parameters needed to define the covariance matrices as in AICC2012 (Bazin et al., CP, 2013 and Veres et al., CP, 2013).
self.lambda_tau=70
self.lambda_a=4000
self.lambda_LID=4000
self.cT1=0.01
self.cT3=0.15
self.sigmam=0.15
