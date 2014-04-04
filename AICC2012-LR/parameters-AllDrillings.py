self.calc_a=False           #Use False for now.
self.calc_tau=False         #Use False for now.
self.calc_LID=False         #Use False for now.
self.calc_a_sigma=False     #Whether sigma on accu is calculated or imported.
self.calc_tau_sigma=False   #Whether sigma on tau is calculated or imported.
self.calc_LID_sigma=False   #Whether sigma on LID is calculated or imported.
self.cT1=0.01
self.cT3=0.15
self.sigmam=0.15
self.corr_tau=np.zeros(51)  #Define the number of nodes of the thinning function.
self.lambda_tau=70
self.lambda_a=4000
self.lambda_LID=4000
self.Dfirn=0.7            #Average density of the firn

