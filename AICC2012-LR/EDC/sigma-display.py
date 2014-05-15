import numpy as np
import matplotlib.pyplot as mpl

readarray=np.loadtxt('accu-sigma-prior.txt')
sigma_a_model_AICC2012=readarray

readarray=np.loadtxt('thinning-sigma-prior.txt')
sigma_tau_model_AICC2012=readarray

readarray=np.loadtxt('LID-sigma-prior.txt')
sigma_LID_model_AICC2012=readarray

readarray=np.loadtxt('output.txt')
depth=readarray[:,0]
sigma_a_model_IC=readarray[:,14]
sigma_tau_model_IC=readarray[:,16]
sigma_LID_model_IC=readarray[:,18]

mpl.figure('Sigma accu')
mpl.plot(sigma_a_model_AICC2012, depth[:-1], label='AICC2012')
mpl.plot(sigma_a_model_IC[:-1], depth[:-1], label='IceChrono')
x1,x2,y1,y2 = mpl.axis()
mpl.axis((x1, x2, y2, 0.))
mpl.legend()

mpl.figure('Sigma thinning')
mpl.plot(sigma_tau_model_AICC2012, depth[:-1], label='AICC2012')
mpl.plot(sigma_tau_model_IC[:-1], depth[:-1], label='IceChrono')
x1,x2,y1,y2 = mpl.axis()
mpl.axis((x1, x2, y2, 0.))
mpl.legend()

mpl.figure('Sigma LID')
mpl.plot(sigma_LID_model_AICC2012, depth, label='AICC2012')
mpl.plot(sigma_LID_model_IC, depth, label='IceChrono')
x1,x2,y1,y2 = mpl.axis()
mpl.axis((x1, x2, y2, 0.))
mpl.legend()

mpl.show()
